#Script for GSE10391 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter


library(GEOquery)
library(limma)

#At OHSU wd
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs")
#My laptop wd
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/GSE10391/")


##########################

#  GSE10391

##########################

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
g.10391.geo <- getGEO("GSE10391", destdir="./") #Downloads the soft files
g.10391.pheno <- phenoData(g.10391.geo[[1]]) #Gets phenotype data
g.10391.pdata <- pData(g.10391.pheno) #Dataframe of all data associated with arrays

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
g.10391.fnames <- as.character(g.10391.pdata$supplementary_file)
g.10391.fnames <- sapply(g.10391.fnames, 
                        function(x) {strsplit(x, "/", fixed=T)[[1]][11]},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
g.10391.arraynames <- as.character(g.10391.pdata$geo_accession)

#Get the treatments for each array
g.10391.treatment <- as.character(g.10391.pdata$description.1)

#Function to parse out the treatment from the description
text_to_label <- function(row){
  #Control logic to separate control from treatments
  if (grepl("Control", row, fixed=T)){
    row.split <- strsplit(row, " ", fixed=T)
    n.days <- row.split[[1]][8]
    label <- paste("Control", n.days, "days", sep="." )
    return(label)
  }
  else {
    row.split <- strsplit(row, " ", fixed=T)
    n.days <- row.split[[1]][17]
    label <- paste("Wayne", n.days, "days", sep=".")
    return(label)
  }
}

#use the filter function to assign treatments
g.10391.treatment <- sapply(g.10391.treatment, text_to_label, USE.NAMES=F)


#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
g.10391.targets <- data.frame(FileName=g.10391.fnames,
                             Cy3 = rep("control", length(g.10391.fnames)),
                             Cy5 = g.10391.treatment)

rownames(g.10391.targets) <- g.10391.arraynames

##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
g.10391.cols <- list(R="Ch2 Intensity (Mean)", G="Ch1 Intensity (Mean)",
                    Rb="Ch2 Background (Median)", Gb="Ch1 Background (Median)")

#Read in the raw intensities. Note that the annotation will be used later
#in the background correction and normalization
#and correspond to Block ,row, and column.
g.10391.rg <- read.maimages(g.10391.targets, columns=g.10391.cols,
                           annotation=c("Sector", 
                                        "X Grid Coordinate (within sector)",
                                        "Y Grid Coordinate (within sector)", 
                                        "Spot", "Name"),
                           path="./GSE10391_RAW")
colnames(g.10391.rg$genes) <- c("Block", "Row", "Column", "Spot", "Name")

#Generate a Layout object for the background correction and normalization
g.10391.rg$printer <- getLayout(g.10391.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./QA/prenormMA")
plotMA3by2(g.10391.rg, path="./QA/prenormMA")

#Normalize the arrays, may have to remove some if the artifacts remain
g.10391.bc <- backgroundCorrect(g.10391.rg, method="normexp", offset=50)
g.10391.bc.norm <- normalizeWithinArrays(g.10391.bc, method="loess")

#Need to filter so that only RV genes are present
g.10391.rv.idx <- grepl(pattern="RV", x=g.10391.rg$genes$Name, fixed=T)
sum(g.10391.rv.idx) #3924 genes represented

g.10391.bc.norm.rv <- g.10391.bc.norm[g.10391.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(g.10391.bc.norm.rv, path="./QA/postnormMA_RVfiltered")

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(g.10391.rg)
dev.off()

png("corrected_densities.png")
plotDensities(g.10391.bc.norm.rv)
dev.off()

################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

#Subset to include only treatment arrays, only genes
g.10391.wayne.idx <- grepl("Wayne", g.10391.targets$Cy5, fixed=T)
g.10391.bc.norm.rv.wayne <- g.10391.bc.norm.rv[, g.10391.wayne.idx]

g.10391.rv.treat.M <- as.data.frame(g.10391.bc.norm.rv.wayne$M)
row.names(g.10391.rv.treat.M) <- g.10391.bc.norm.rv.wayne$genes$Name

discretizer <- function(val){
  if(val >= 1){
    return(1)
  }
  
  if(val <= -1){
    return (-1)
  }
  return(0)
}

vec.discret <- function(vec){
  return(sapply(vec, discretizer, USE.NAMES=F))
}

g.10391.disc <- data.frame(lapply(g.10391.rv.treat.M, vec.discret),
                          row.names = rownames(g.10391.rv.treat.M))

save(g.10391.disc, file="../GSE10391_disc.RData")
save(g.10391.rv.treat.M, file="../GSE10391_M.RData")
