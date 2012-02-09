#Script for GSE16146 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter

##TODO
##8523 has the reference as Cy5 for some samples

library(GEOquery)
library(limma)

#At OHSU Dropbox
setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
#setwd("~/schoolDB/Dropbox/thesis_work")
source("./code/exprs/exprs_funcs.R")

setwd("./data/exprs/GSE16146")

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
g.16146.geo <- getGEO("GSE16146", destdir="./", GSEMatrix=F) #Downloads the soft files
g.16146.geo.matrix <- getGEO("GSE16146", destdir="./")

#Need to pull out the data table for each array and save as .txt for use in read.maimages
#Get a list of the array objects
g.16146.gsm <- GSMList(g.16146.geo)

#load platform data
gpl.8523 <- g.16146.geo.matrix[["GSE16146-GPL8523_series_matrix.txt.gz"]]
gpl.8523.pheno <- phenoData(gpl.8523) #Gets phenotype data
gpl.8523.pdata <- pData(gpl.8523.pheno) #Dataframe of all data associated with arrays

gpl.8523.annot <- featureData(gpl.8523)@data

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

gpl.8523.rg <- read.maimages(gpl.8523.targets,
                             annotation=c("ID_REF"),
                             columns=gpl.8523.cols,
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
plotMA3by2(gpl.8523.rg, path="./QA/prenormMA", main=gpl.8523.targets$Cy5)


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.8523.bc <- backgroundCorrect(gpl.8523.rg, method="normexp", offset=50)
gpl.8523.bc.norm <- normalizeWithinArrays(gpl.8523.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.8523.rv.idx <- grepl(pattern="RV", x=gpl.8523.rg$genes$ORF, fixed=T)
sum(gpl.8523.rv.idx) #16,068 features represented
gpl.8523.bc.norm.rv <- gpl.8523.bc.norm[gpl.8523.rv.idx,]

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.8523.bc.norm.rv, path="./QA/postnormMA_RVfiltered", main=gpl.8523.targets$Cy5)

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

gpl.8523.rv.M <- as.data.frame(gpl.8523.bc.norm.rv$M)
gpl.8523.disc <- discretize(gpl.8523.rv.M)
gpl.8523.disc$gene <- gpl.8523.bc.norm.rv$genes$ORF


#######
# 
# GPL 4293
#
#######

gpl.4293 <- g.16146.geo[["GSE16146-GPL4293_series_matrix.txt.gz"]]
gpl.4293.pheno <- phenoData(gpl.4293) #Gets phenotype data
gpl.4293.pdata <- pData(gpl.4293.pheno) #Dataframe of all data associated with arrays

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.4293.fnames <- as.character(gpl.4293.pdata$supplementary_file)
gpl.4293.fnames <- sapply(gpl.4293.fnames, 
                        function(x) {strsplit(x, "/", fixed=T)[[1]][11]},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.4293.arraynames <- as.character(gpl.4293.pdata$geo_accession)

#Get the treatments for each array
gpl.4293.treatment <- as.character(gpl.4293.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.4293.targets <- data.frame(FileName=gpl.4293.fnames,
                             Cy3 = rep("control", length(gpl.4293.fnames)),
                             Cy5 = gpl.4293.treatment, stringsAsFactors=F)

rownames(gpl.4293.targets) <- gpl.4293.arraynames

#Initial attempts at reading in the images threw errors on multiple arrays.
#Inspection of the files shows they are corrupted in some way.
#Exclude from targets dataframe.
dim(gpl.4293.targets) #12 x 3
bad.arrays <- c("GSM237676.gpr.gz")
bad.idx <- which(gpl.4293.targets$FileName %in% bad.arrays)
gpl.4293.targets <- gpl.4293.targets[-bad.idx,]
dim(gpl.4293.targets) #1x3, 1 arrays removed successfully


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.4293.cols <- list(R="F635 Mean", G="F532 Mean",
                    Rb="B635 Median", Gb="B532 Median")

gpl.4293.rg <- read.maimages(gpl.4293.targets,
                             source="genepix",
                             columns=gpl.4293.cols,
                             path="./GSE16146_RAW")

#Generate a Layout object for the background correction and normalization
gpl.4293.rg$printer <- getLayout(gpl.4293.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL4293/QA/prenormMA", recursive=T)
plotMA3by2(gpl.4293.rg, path="./GPL4293/QA/prenormMA")


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.4293.bc <- backgroundCorrect(gpl.4293.rg, method="normexp", offset=50)
gpl.4293.bc.norm <- normalizeWithinArrays(gpl.4293.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.4293.rv.idx <- grepl(pattern="Rv", x=gpl.4293.rg$genes$Name, fixed=T)
sum(gpl.4293.rv.idx) #24,126 probes
gpl.4293.bc.norm.rv <- gpl.4293.bc.norm[gpl.4293.rv.idx,]
length(unique(gpl.4293.bc.norm.rv$genes$Name)) #3900 genes

#Realizing that there are multiple spots for each gene. Will need to average
#the ratios or deal with it in some way.

#Redo the MA plots and see if artifacts disappear
dir.create("./GPL4293/QA/postnormMA_RVfiltered")
plotMA3by2(gpl.4293.bc.norm.rv, path="./GPL4293/QA/postnormMA_RVfiltered")

#QA plots to see if normalization worked
setwd("./GPL4293/QA")
png("uncorrected_densities.png")
plotDensities(gpl.4293.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.4293.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

gpl.4293.rv.M <- as.data.frame(gpl.4293.bc.norm.rv$M)
gpl.4293.rv.M$gene <- gpl.4293.bc.norm.rv$genes$Name
gpl.4293.gene_ids <- unique(gpl.4293.rv.M$gene)
gpl.4293.M.avg <- avg_probes(gpl.4293.rv.M, gpl.4293.gene_ids)

gpl.4293.disc <- discretize(gpl.4293.M.avg)


#######
# 
# GPL 5774
#
#######

gpl.5774 <- g.16146.geo[["GSE16146-GPL5774_series_matrix.txt.gz"]]
gpl.5774.pheno <- phenoData(gpl.5774) #Gets phenotype data
gpl.5774.pdata <- pData(gpl.5774.pheno) #Dataframe of all data associated with arrays

#Pull the file names out of the pData dataframe, need to process the strings a
#bit because they point to the ftp path to download, we only care about file name
gpl.5774.fnames <- as.character(gpl.5774.pdata$supplementary_file)
gpl.5774.fnames <- sapply(gpl.5774.fnames, 
                        function(x) {strsplit(x, "/", fixed=T)[[1]][11]},
                        USE.NAMES=F
                        )

#Get the GEO accessions for each array to use as an ID
gpl.5774.arraynames <- as.character(gpl.5774.pdata$geo_accession)

#Get the treatments for each array
gpl.5774.treatment <- as.character(gpl.5774.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
gpl.5774.targets <- data.frame(FileName=gpl.5774.fnames,
                             Cy3 = rep("control", length(gpl.5774.fnames)),
                             Cy5 = gpl.5774.treatment, stringsAsFactors=F)

rownames(gpl.5774.targets) <- gpl.5774.arraynames

#Initial attempts at reading in the images threw errors on multiple arrays.
#Inspection of the files shows they are corrupted in some way.
#Exclude from targets dataframe.
dim(gpl.5774.targets) #12 x 3
bad.arrays <- c("GSM237676.gpr.gz")
bad.idx <- which(gpl.5774.targets$FileName %in% bad.arrays)
gpl.5774.targets <- gpl.5774.targets[-bad.idx,]
dim(gpl.5774.targets) #1x3, 1 arrays removed successfully


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.5774.cols <- list(R="F635 Mean", G="F532 Mean",
                    Rb="B635 Median", Gb="B532 Median")

gpl.5774.rg <- read.maimages(gpl.5774.targets,
                             source="genepix",
                             columns=gpl.5774.cols,
                             annotation=c("Block","Row","Column","ID","NAME"),
                             path="./GSE16146_RAW")

#Generate a Layout object for the background correction and normalization
gpl.5774.rg$printer <- getLayout(gpl.5774.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL5774/QA/prenormMA", recursive=T)
plotMA3by2(gpl.5774.rg, path="./GPL5774/QA/prenormMA")


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.5774.bc <- backgroundCorrect(gpl.5774.rg, method="normexp", offset=50)
gpl.5774.bc.norm <- normalizeWithinArrays(gpl.5774.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.5774.rv.idx <- grepl(pattern="Rv", x=gpl.5774.rg$genes$NAME, fixed=T)
sum(gpl.5774.rv.idx) #16,076 probes
gpl.5774.bc.norm.rv <- gpl.5774.bc.norm[gpl.5774.rv.idx,]
length(unique(gpl.5774.bc.norm.rv$genes$NAME)) #3898 genes

#Redo the MA plots and see if artifacts disappear
dir.create("./GPL5774/QA/postnormMA_RVfiltered")
plotMA3by2(gpl.5774.bc.norm.rv, path="./GPL5774/QA/postnormMA_RVfiltered")

#QA plots to see if normalization worked
setwd("./GPL5774/QA")
png("uncorrected_densities.png")
plotDensities(gpl.5774.rg)
dev.off()

png("corrected_densities.png")
plotDensities(gpl.5774.bc.norm.rv)
dev.off()


################
#
# Begin discretizing values
#
################

##Extract log-2 expression ratios and discretize

gpl.5774.rv.M <- as.data.frame(gpl.5774.bc.norm.rv$M)
gpl.5774.rv.M$gene <- gpl.5774.bc.norm.rv$genes$NAME
gpl.5774.gene_ids <- unique(gpl.5774.rv.M$gene)
gpl.5774.M.avg <- avg_probes(gpl.5774.rv.M, gpl.5774.gene_ids)

gpl.5774.disc <- discretize(gpl.5774.M.avg)

