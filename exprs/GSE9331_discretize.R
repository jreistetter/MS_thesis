#Script for GSE9331 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter

library(GEOquery)
library(limma)

#At OHSU wd
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs/GSE9331")
#My laptop wd
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/GSE9331/")

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
g.9331.geo <- getGEO("GSE9331", destdir="./") #Downloads the soft files

#This dataset has 3 different platforms, so work with each separately

#######
# 
# GPL 4291
#
#######

gpl.4291 <- g.9331.geo[["GSE9331-GPL4291_series_matrix.txt.gz"]]
gpl.4291.pheno <- phenoData(gpl.4291) #Gets phenotype data
gpl.4291.pdata <- pData(gpl.4291.pheno) #Dataframe of all data associated with arrays

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.4291.fnames <- as.character(gpl.4291.pdata$supplementary_file)
gpl.4291.fnames <- sapply(gpl.4291.fnames, 
                        function(x) {strsplit(x, "/", fixed=T)[[1]][11]},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.4291.arraynames <- as.character(gpl.4291.pdata$geo_accession)

#Get the treatments for each array
gpl.4291.treatment <- as.character(gpl.4291.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.4291.targets <- data.frame(FileName=gpl.4291.fnames,
                             Cy3 = rep("control", length(gpl.4291.fnames)),
                             Cy5 = gpl.4291.treatment, stringsAsFactors=F)

rownames(gpl.4291.targets) <- gpl.4291.arraynames

#Initial attempts at reading in the images threw errors on multiple arrays.
#Inspection of the files shows they are corrupted in some way.
#Exclude from targets dataframe.
dim(gpl.4291.targets) #38 x 3
bad.arrays <- c("GSM237638.gpr.gz", "GSM237639.gpr.gz", 
                "GSM237640.gpr.gz", "GSM237641.gpr.gz",
                "GSM237642.gpr.gz", "GSM237647.gpr.gz",
                "GSM237648.gpr.gz")
bad.idx <- which(gpl.4291.targets$FileName %in% bad.arrays)
gpl.4291.targets <- gpl.4291.targets[-bad.idx,]
dim(gpl.4291.targets) #31x3, 7 arrays removed successfully


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.4291.cols <- list(R="F635 Mean", G="F532 Mean",
                    Rb="B635 Median", Gb="B532 Median")

gpl.4291.rg <- read.maimages(gpl.4291.targets,
                             source="genepix",
                             columns=gpl.4291.cols,
                             path="./GSE9331_RAW")

#Generate a Layout object for the background correction and normalization
gpl.4291.rg$printer <- getLayout(gpl.4291.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./QA")
dir.create("./QA/prenormMA")
plotMA3by2(gpl.4291.rg, path="./QA/prenormMA")


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.4291.bc <- backgroundCorrect(gpl.4291.rg, method="normexp", offset=50)
gpl.4291.bc.norm <- normalizeWithinArrays(gpl.4291.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.4291.rv.idx <- grepl(pattern="Rv", x=gpl.4291.rg$genes$Name, fixed=T)
sum(gpl.4291.rv.idx) #16068 genes represented
gpl.4291.bc.norm.rv <- gpl.4291.bc.norm[gpl.4291.rv.idx,]

#Realizing that there are multiple spots for each gene. Will need to average
#the ratios or deal with it in some way.

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.4291.bc.norm.rv, path="./QA/postnormMA_RVfiltered")

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(gpl.4291.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.4291.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

gpl.4291.rv.M <- as.data.frame(gpl.4291.bc.norm.rv$M)
row.names(gpl.4291.rv.M) <- gpl.4291.bc.norm.rv$genes$Name

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

gpl.4291.disc <- data.frame(lapply(gpl.4291.rv.M, vec.discret))



