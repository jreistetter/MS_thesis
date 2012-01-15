#Script to download the supplementary files (raw data) for 5 GEO datasets 
#that will be used to build the PMN.

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

g.8786.rg <- read.maimages(g.8786.targets, columns=g.8786.cols,
                           annotation=c("Sector", "X Grid Coordinate (within sector)",
                                        "Y Grid Coordinate (within sector)", 
                                        "Spot", "Name"),
                           path="./GSE8786_RAW")
colnames(g.8786.rg$genes) <- c("Block", "Row", "Column", "Spot", "Name")

g.8786.rg$printer <- getLayout(g.8786.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./QA/prenormMA")
plotMA3by2(g.8786.rg, path="./QA/prenormMA")



#The filtered helped remove some of the artifacts, but some odd
#patterns (diagonal lines fanning out from the lower left quadrant
#of the plots)

#Normalize the arrays, may have to remove some if the artifacts remain
g.8786.bc <- backgroundCorrect(g.8786.rg, method="normexp", offset=50)
g.8786.bc.norm <- normalizeWithinArrays(g.8786.bc)
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

plotDensities(g.8786.MAq.rv)

################
#
# Begin differential expression analysis
#
################

#all arrays
fit <- lmFit(g.8786.bc.norm.rv)
fit <- eBayes(fit)
topTable(fit)

#just wayne growth arrays
g.8786.wayne.idx <- grepl("Wayne", g.8786.targets$Cy5, fixed=T)
g.8786.bc.norm.rv.wayne <- g.8786.bc.norm.rv[, g.8786.wayne.idx]

fit.wayne <- lmFit(g.8786.bc.norm.rv.wayne)
fit.wayne <- eBayes(fit.wayne)
topTable(fit.wayne)

