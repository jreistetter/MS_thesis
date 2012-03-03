#Script for GSE16146 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter


library(GEOquery)
library(limma)

#At OHSU Dropbox
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
setwd("~/schoolDB/Dropbox/thesis_work")
source("./code/exprs/exprs_funcs.R")

setwd("./data/exprs/GSE16146")

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
g.16146.geo <- getGEO("GSE16146", destdir="./", GSEMatrix=F) #Downloads the soft files
g.16146.geo.matrix <- getGEO("GSE16146", destdir="./")

#Need to pull out the data table for each array and save as .txt for use in read.maimages
#Get a list of the array objects
g.16146.gsm <- GSMList(g.16146.geo)

#Function to pull out the data table and save as text
save_table <- function(gsm){
  dat <- Table(gsm@dataTable)
  accession <- gsm@header$geo_accession
  fname <- paste(accession, ".txt", sep="")
  write.table(dat, file=fname, sep="\t", quote=F, col.names=T, row.names=F)
  
}

dir.create("./GSE16146_RAW")
setwd("./GSE16146_RAW")

#Save out to text
lapply(g.16146.gsm, save_table)


#This dataset has 3 different platforms, so work with each separately
setwd("..")

#######
# 
# GPL 8523
#
#######
#load platform data
gpl.8523 <- g.16146.geo.matrix[["GSE16146-GPL8523_series_matrix.txt.gz"]]
gpl.8523.pheno <- phenoData(gpl.8523) #Gets phenotype data
gpl.8523.pdata <- pData(gpl.8523.pheno) #Dataframe of all data associated with arrays

gpl.8523.annot <- featureData(gpl.8523)@data

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.8523.fnames <- rownames(gpl.8523.pdata)
gpl.8523.fnames <- sapply(gpl.8523.fnames, 
                        function(x) {paste(x, ".txt", sep="")},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.8523.arraynames <- rownames(gpl.8523.pdata)

#Get the treatments for each array
gpl.8523.treatment <- as.character(gpl.8523.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.8523.targets <- data.frame(FileName=gpl.8523.fnames,
                             Cy3 = rep("control", length(gpl.8523.fnames)),
                             Cy5 = gpl.8523.treatment, stringsAsFactors=F)

rownames(gpl.8523.targets) <- gpl.8523.arraynames


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.8523.cols <- list(R="CH2I_MEAN", G="CH1I_MEAN",
                    Rb="CH2B_MEDIAN", Gb="CH1B_MEDIAN")

#Define weight function to exclude spots that have a "FLAG" column value < 0 in the
#raw data
flagged <- function(x) as.numeric(x[,"FLAG"] > -1)

gpl.8523.rg <- read.maimages(gpl.8523.targets,
                             annotation=c("ID_REF"),
                             columns=gpl.8523.cols,
                             wt.fun=flagged,
                             path="./GSE16146_RAW")

gpl.8523.gal <- readGAL("SMD_print_745.gal")
colnames(gpl.8523.gal)[5] <- "Reporter.Name"
#Generate a Layout object for the background correction and normalization
gpl.8523.rg$genes <- merge(gpl.8523.annot, gpl.8523.rg$genes, by.x = "ID", by.y="ID_REF", all.y=T)
#gpl.8523.rg$printer <- getLayout(gpl.8523.rg$genes)
gpl.8523.rg$printer <- getLayout(gpl.8523.gal)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL8523")
setwd("./GPL8523")
dir.create("./QA")
dir.create("./QA/prenormMA")
plotMA3by2(gpl.8523.rg, path="./QA/prenormMA", 
           main=paste(gpl.8523.targets$FileName, gpl.8523.targets$Cy5, sep=" - "))


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.8523.bc <- backgroundCorrect(gpl.8523.rg, method="normexp", offset=50)
gpl.8523.bc.norm <- normalizeWithinArrays(gpl.8523.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.8523.rv.idx <- grepl(pattern="RV", x=gpl.8523.rg$genes$ORF, fixed=T)
sum(gpl.8523.rv.idx) #16,068 features represented
gpl.8523.bc.norm.rv <- gpl.8523.bc.norm[gpl.8523.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.8523.bc.norm.rv, path="./QA/postnormMA_RVfiltered", 
           main=paste(gpl.8523.targets$FileName, gpl.8523.targets$Cy5, sep=" - "))

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(gpl.8523.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.8523.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize
#override remove_bad_spots function because gene name is
#in different named column for GSE16146
remove_bad_spots <- function(ma_list){
  probe.weights <- ma_list$weights
  probe.weights[probe.weights == 0] <- NA
  cleaned <- probe.weights * ma_list$M
  rownames(cleaned) <- ma_list$genes$ORF
  return(as.data.frame(cleaned))
  
}
gpl.8523.rv.M <- remove_bad_spots(gpl.8523.bc.norm.rv)
gpl.8523.disc <- discretize(gpl.8523.rv.M)


#######
# 
# GPL 8561
#
#######
setwd("../..")

#load platform data
gpl.8561 <- g.16146.geo.matrix[["GSE16146-GPL8561_series_matrix.txt.gz"]]
gpl.8561.pheno <- phenoData(gpl.8561) #Gets phenotype data
gpl.8561.pdata <- pData(gpl.8561.pheno) #Dataframe of all data associated with arrays

gpl.8561.annot <- featureData(gpl.8561)@data

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.8561.fnames <- rownames(gpl.8561.pdata)
gpl.8561.fnames <- sapply(gpl.8561.fnames, 
                        function(x) {paste(x, ".txt", sep="")},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.8561.arraynames <- rownames(gpl.8561.pdata)

#Get the treatments for each array
gpl.8561.treatment <- as.character(gpl.8561.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.8561.targets <- data.frame(FileName=gpl.8561.fnames,
                             Cy3 = rep("control", length(gpl.8561.fnames)),
                             Cy5 = gpl.8561.treatment, stringsAsFactors=F)

rownames(gpl.8561.targets) <- gpl.8561.arraynames


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.8561.cols <- list(R="CH2I_MEAN", G="CH1I_MEAN",
                    Rb="CH2B_MEDIAN", Gb="CH1B_MEDIAN")

gpl.8561.rg <- read.maimages(gpl.8561.targets,
                             annotation=c("ID_REF"),
                             columns=gpl.8561.cols,
                             wt.fun=flagged,
                             path="./GSE16146_RAW")

gpl.8561.gal <- readGAL("SMD_print_498.gal")
colnames(gpl.8561.gal)[5] <- "Reporter.Name"
#Generate a Layout object for the background correction and normalization
gpl.8561.rg$genes <- merge(gpl.8561.annot, gpl.8561.rg$genes, by.x = "ID", by.y="ID_REF", all.y=T)
#gpl.8561.rg$printer <- getLayout(gpl.8561.rg$genes)
gpl.8561.rg$printer <- getLayout(gpl.8561.gal)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL8561")
setwd("./GPL8561")
dir.create("./QA")
dir.create("./QA/prenormMA")
plotMA3by2(gpl.8561.rg, path="./QA/prenormMA", 
           main=paste(gpl.8561.targets$FileName, gpl.8561.targets$Cy5, paste=" - "))


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.8561.bc <- backgroundCorrect(gpl.8561.rg, method="normexp", offset=50)
gpl.8561.bc.norm <- normalizeWithinArrays(gpl.8561.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.8561.rv.idx <- grepl(pattern="RV", x=gpl.8561.rg$genes$ORF, fixed=T)
sum(gpl.8561.rv.idx) #4,686 features represented
gpl.8561.bc.norm.rv <- gpl.8561.bc.norm[gpl.8561.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.8561.bc.norm.rv, path="./QA/postnormMA_RVfiltered", 
           main=paste(gpl.8561.targets$FileName, gpl.8561.targets$Cy5, paste=" - "))

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(gpl.8561.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.8561.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

gpl.8561.rv.M <- remove_bad_spots(gpl.8561.bc.norm.rv)
gpl.8561.disc <- discretize(gpl.8561.rv.M.avg)


#######
# 
# GPL 8562
#
#######
setwd("../..")

#load platform data
gpl.8562 <- g.16146.geo.matrix[["GSE16146-GPL8562_series_matrix.txt.gz"]]
gpl.8562.pheno <- phenoData(gpl.8562) #Gets phenotype data
gpl.8562.pdata <- pData(gpl.8562.pheno) #Dataframe of all data associated with arrays

gpl.8562.annot <- featureData(gpl.8562)@data

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.8562.fnames <- rownames(gpl.8562.pdata)
gpl.8562.fnames <- sapply(gpl.8562.fnames, 
                        function(x) {paste(x, ".txt", sep="")},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.8562.arraynames <- rownames(gpl.8562.pdata)

#Get the treatments for each array
gpl.8562.treatment <- as.character(gpl.8562.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.8562.targets <- data.frame(FileName=gpl.8562.fnames,
                             Cy3 = rep("control", length(gpl.8562.fnames)),
                             Cy5 = gpl.8562.treatment, stringsAsFactors=F)

rownames(gpl.8562.targets) <- gpl.8562.arraynames


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.8562.cols <- list(R="CH2I_MEAN", G="CH1I_MEAN",
                    Rb="CH2B_MEDIAN", Gb="CH1B_MEDIAN")

gpl.8562.rg <- read.maimages(gpl.8562.targets,
                             annotation=c("ID_REF"),
                             columns=gpl.8562.cols,
                             wt.fun=flagged,
                             path="./GSE16146_RAW")

gpl.8562.gal <- readGAL("SMD_print_580.gal")
colnames(gpl.8562.gal)[5] <- "Reporter.Name"
#Generate a Layout object for the background correction and normalization
gpl.8562.rg$genes <- merge(gpl.8562.annot, gpl.8562.rg$genes, by.x = "ID", by.y="ID_REF", all.y=T)
#gpl.8562.rg$printer <- getLayout(gpl.8562.rg$genes)
gpl.8562.rg$printer <- getLayout(gpl.8562.gal)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL8562")
setwd("./GPL8562")
dir.create("./QA")
dir.create("./QA/prenormMA")
plotMA3by2(gpl.8562.rg, path="./QA/prenormMA", 
           main=paste(gpl.8562.targets$FileName, gpl.8562.targets$Cy5, sep=" - "))


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.8562.bc <- backgroundCorrect(gpl.8562.rg, method="normexp", offset=50)
gpl.8562.bc.norm <- normalizeWithinArrays(gpl.8562.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.8562.rv.idx <- grepl(pattern="RV", x=gpl.8562.rg$genes$ORF, fixed=T)
sum(gpl.8562.rv.idx) #16,068 features represented
gpl.8562.bc.norm.rv <- gpl.8562.bc.norm[gpl.8562.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.8562.bc.norm.rv, path="./QA/postnormMA_RVfiltered", 
           main=paste(gpl.8562.targets$FileName, gpl.8562.targets$Cy5, sep=" - "))

#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(gpl.8562.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.8562.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize
gpl.8562.rv.M <- remove_bad_spots(gpl.8562.bc.norm.rv)
gpl.8562.disc <- discretize(gpl.8562.rv.M)

#Save all the M objects for later
setwd("../..")

save(gpl.8523.rv.M, file="gpl.8523.M.RData")
save(gpl.8561.rv.M, file="gpl.8561.M.RData")
save(gpl.8562.rv.M, file="gpl.8562.M.RData")



