#Script for GSE8786 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter

#Notes about the datasets:
#GSE8786 - Raw data good, text format
#GSE8839 - appears to be corrupted (some columns shifted in raw data)
            #*Emailing GEO, ignore this dataset for now
#GSE9331 - Raw data good, genepix file format (!) (.gpr)
          #It's 3 different platforms, and 4921 says channel 1 is cy5.
          #Need to check this later to see if just entry error in GEO
          #gpl4291 ch1 = cy5 = red
          #gpl4293 ch1 = cy3 = green
          #gpl5774 ch1 = cy3 = green
#GSE10391 - Raw data good, text format
#GSE16146 - get softfile via getGEO, raw data not available separately


library(GEOquery)
library(limma)

#At OHSU wd
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs")
#My laptop wd
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/GSE8786/")


##########################

#  GSE8786

##########################

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
g.8786.geo <- getGEO("GSE8786", destdir="./") #Downloads the soft files
g.8786.pheno <- phenoData(g.8786.geo[[1]]) #Gets phenotype data
g.8786.pdata <- pData(g.8786.pheno) #Dataframe of all data associated with arrays

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
g.8786.fnames <- as.character(g.8786.pdata$supplementary_file)
g.8786.fnames <- sapply(g.8786.fnames, 
                        function(x) {strsplit(x, "/", fixed=T)[[1]][11]},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
g.8786.arraynames <- as.character(g.8786.pdata$geo_accession)

#Get the treatments for each array
g.8786.treatment <- as.character(g.8786.pdata$description.1)

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
g.8786.treatment <- sapply(g.8786.treatment, text_to_label, USE.NAMES=F)


#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
g.8786.targets <- data.frame(FileName=g.8786.fnames,
                             Cy3 = rep("control", length(g.8786.fnames)),
                             Cy5 = g.8786.treatment)

rownames(g.8786.targets) <- g.8786.arraynames

##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
g.8786.cols <- list(R="Ch2 Intensity (Mean)", G="Ch1 Intensity (Mean)",
                    Rb="Ch2 Background (Median)", Gb="Ch1 Background (Median)")

#Read in the raw intensities. Note that the annotation will be used later
#in the background correction and normalization
#and correspond to Block ,row, and column.
g.8786.rg <- read.maimages(g.8786.targets, columns=g.8786.cols,
                           annotation=c("Sector", 
                                        "X Grid Coordinate (within sector)",
                                        "Y Grid Coordinate (within sector)", 
                                        "Spot", "Name"),
                           path="./GSE8786_RAW")
colnames(g.8786.rg$genes) <- c("Block", "Row", "Column", "Spot", "Name")

#Generate a Layout object for the background correction and normalization
g.8786.rg$printer <- getLayout(g.8786.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./QA/prenormMA")
plotMA3by2(g.8786.rg, path="./QA/prenormMA")

#Normalize the arrays, may have to remove some if the artifacts remain
g.8786.bc <- backgroundCorrect(g.8786.rg, method="normexp", offset=50)
g.8786.bc.norm <- normalizeWithinArrays(g.8786.bc, method="loess")

#Need to filter so that only RV genes are present
g.8786.rv.idx <- grepl(pattern="RV", x=g.8786.rg$genes$Name, fixed=T)
sum(g.8786.rv.idx) #3924 genes represented

g.8786.bc.norm.rv <- g.8786.bc.norm[g.8786.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./postnormMA_RVfiltered")
plotMA3by2(g.8786.bc.norm.rv, path="./postnormMA_RVfiltered")

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(g.8786.rg)
dev.off()

png("corrected_densities.png")
plotDensities(g.8786.bc.norm.rv)
dev.off()

################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

#Subset to include only treatment arrays, only genes
g.8786.wayne.idx <- grepl("Wayne", g.8786.targets$Cy5, fixed=T)
g.8786.bc.norm.rv.wayne <- g.8786.bc.norm.rv[, g.8786.wayne.idx]

g.8786.rv.treat.M <- as.data.frame(g.8786.bc.norm.rv.wayne$M)
row.names(g.8786.rv.treat.M) <- g.8786.bc.norm.rv.wayne$genes$Name

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

g.8786.disc <- data.frame(lapply(g.8786.rv.treat.M, vec.discret),
                          row.names = rownames(g.8786.rv.treat.M))
