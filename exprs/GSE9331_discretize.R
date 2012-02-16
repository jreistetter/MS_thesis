#Script for GSE9331 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter

##TODO
##4291 has the reference as Cy5 for some samples

library(GEOquery)
library(limma)

#At OHSU Dropbox
setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
#setwd("~/schoolDB/Dropbox/thesis_work")
source("./code/exprs/exprs_funcs.R")

setwd("./data/exprs/GSE9331")

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
                             wt.fun=wtflags(weight=0, cutoff=-1),
                             path="./GSE9331_RAW")

#Generate a Layout object for the background correction and normalization
gpl.4291.rg$printer <- getLayout(gpl.4291.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL4291")
setwd("./GPL4291")
dir.create("./QA")
dir.create("./QA/prenormMA")
plotMA3by2(gpl.4291.rg, path="./QA/prenormMA", 
           main=paste(gpl.4291.targets$FileName, gpl.4291.targets$Cy5, sep=" - "))


#Normalize the arrays, may have to remove some if the artifacts remain
gpl.4291.bc <- backgroundCorrect(gpl.4291.rg, method="normexp", offset=50)
gpl.4291.bc.norm <- normalizeWithinArrays(gpl.4291.bc, method="loess")

#Need to filter so that only RV genes are present
gpl.4291.rv.idx <- grepl(pattern="Rv", x=gpl.4291.rg$genes$Name, fixed=T)
sum(gpl.4291.rv.idx) #16,068 features represented
gpl.4291.bc.norm.rv <- gpl.4291.bc.norm[gpl.4291.rv.idx,]
length(unique(gpl.4291.bc.norm.rv$genes$Name)) #4595 genes

#Redo the MA plots and see if artifacts disappear
dir.create("./QA/postnormMA_RVfiltered")
plotMA3by2(gpl.4291.bc.norm.rv, path="./QA/postnormMA_RVfiltered",
           main=paste(gpl.4291.targets$FileName, gpl.4291.targets$Cy5, sep=" - "))

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
gpl.4291.rv.M <- remove_bad_spots(gpl.4291.bc.norm.rv)
gpl.4291.rv.M$gene <- rownames(gpl.4291.rv.M)
gpl.4291.gene_ids <- unique(gpl.4291.rv.M$gene)

gpl.4291.M.avg <- avg_probes(gpl.4291.rv.M, gpl.4291.gene_ids)

##Need to reverse the sign of arrays where the control was Cy5
ch2_Cy3.idx <- which(gpl.4291.pdata$label_ch2 == "Cy3")
gpl.4291.ch2_Cy3.arrays <- as.character(gpl.4291.pdata$geo_accession[ch2_Cy3.idx])
for (i in 1:length(gpl.4291.ch2_Cy3.arrays)){
  gpl.4291.ch2_Cy3.arrays[i] <- paste(gpl.4291.ch2_Cy3.arrays[i], ".gpr", sep="")
}
gpl.4291.ch2_Cy3.idx <- which(colnames(gpl.4291.M.avg) %in% gpl.4291.ch2_Cy3.arrays)

for (idx in gpl.4291.ch2_Cy3.idx){
  gpl.4291.M.avg[,idx] <- -1 * gpl.4291.M.avg[,idx]
}


gpl.4291.disc <- discretize(gpl.4291.M.avg)

#Get coefficient of variation of ratios of the the replicated probes
#Remove bad arrays first
bad.4291 <- c("GSM237637.gpr", "GSM237644.gpr", "GSM237645.gpr", "GSM237646.gpr", "GSM237650.gpr", 
              "GSM237651.gpr", "GSM237655.gpr", "GSM237666.gpr", "GSM237667.gpr")

good.4291.col.idx <- which(!(colnames(gpl.4291.rv.M) %in% bad.4291))
gpl.4291.clean <- gpl.4291.rv.M[,good.4291.col.idx]

replicated.4291 <- unique(gpl.4291.clean[duplicated(gpl.4291.clean$gene),]$gene)
coefs.4291 <- sapply(replicated.4291, probe_CV, df=gpl.4291.rv.M)
coefs.mean.4291 <- apply(coefs.4291, 2, mean)
coefs.median.4291 <- apply(coefs.4291, 2, median)

png("GPL4291_probe_mean_CoV.png")
boxplot(coefs.mean.4291, 
        main = "Gene-wise mean coefficient of variation\nProbe log(R/G)\nGPL4291 22 arrays 3897 genes",
        ylab="Mean CoV")
dev.off()

png("GPL4291_probe_median_CoV.png")
boxplot(coefs.median.4291, 
        main = "Gene-wise median coefficient of variation\nProbe log(R/G)\nGPL4291 22 arrays 3897 genes",
        ylab="Median CoV")
dev.off()

#######
# 
# GPL 4293
#
#######
setwd("../../")

gpl.4293 <- g.9331.geo[["GSE9331-GPL4293_series_matrix.txt.gz"]]
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
                             wt.fun=wtflags(weight=0, cutoff=-1),
                             path="./GSE9331_RAW")

#Generate a Layout object for the background correction and normalization
gpl.4293.rg$printer <- getLayout(gpl.4293.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL4293/QA/prenormMA", recursive=T)
plotMA3by2(gpl.4293.rg, path="./GPL4293/QA/prenormMA", 
           main=paste(gpl.4293.targets$FileName, gpl.4293.targets$Cy5, sep=" - "))


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
plotMA3by2(gpl.4293.bc.norm.rv, path="./GPL4293/QA/postnormMA_RVfiltered", 
           main=paste(gpl.4293.targets$FileName, gpl.4293.targets$Cy5, sep=" - "))

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
gpl.4293.rv.M <- remove_bad_spots(gpl.4293.bc.norm.rv)
gpl.4293.rv.M$gene <- rownames(gpl.4293.rv.M)
gpl.4293.gene_ids <- unique(gpl.4293.rv.M$gene)

gpl.4293.M.avg <- avg_probes(gpl.4293.rv.M, gpl.4293.gene_ids)


replicated.4293 <- unique(gpl.4293.rv.M[duplicated(gpl.4293.rv.M$gene),]$gene)
coefs.4293 <- sapply(replicated.4293, probe_CV, df=gpl.4293.rv.M)
coefs.mean.4293 <- apply(coefs.4293, 2, mean)
coefs.median.4293 <- apply(coefs.4293, 2, median)

png("GPL4293_probe_mean_CoV.png")
boxplot(coefs.mean.4293, 
        main = "Gene-wise mean coefficient of variation\nProbe log(R/G)\nGPL4293 11 arrays 3897 genes",
        ylab="Mean CoV")
dev.off()

png("GPL4293_probe_median_CoV.png")
boxplot(coefs.median.4293, 
        main = "Gene-wise median coefficient of variation\nProbe log(R/G)\nGPL4293 11 arrays 3897 genes",
        ylab="Median CoV")
dev.off()


#######
# 
# GPL 5774
#
#######

setwd("../..")

gpl.5774 <- g.9331.geo[["GSE9331-GPL5774_series_matrix.txt.gz"]]
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


##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
gpl.5774.cols <- list(R="F635 Mean", G="F532 Mean",
                    Rb="B635 Median", Gb="B532 Median")

gpl.5774.rg <- read.maimages(gpl.5774.targets,
                             source="genepix",
                             columns=gpl.5774.cols,
                             annotation=c("Block","Row","Column","ID","NAME"),
                             path="./GSE9331_RAW")

#Generate a Layout object for the background correction and normalization
gpl.5774.rg$printer <- getLayout(gpl.5774.rg$genes)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./GPL5774/QA/prenormMA", recursive=T)
plotMA3by2(gpl.5774.rg, path="./GPL5774/QA/prenormMA", 
           main=paste(gpl.5774.targets$FileName, gpl.5774.targets$Cy5, sep=" - "))


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
plotMA3by2(gpl.5774.bc.norm.rv, path="./GPL5774/QA/postnormMA_RVfiltered", 
           main=paste(gpl.5774.targets$FileName, gpl.5774.targets$Cy5, sep=" - "))

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

#Save M objects for later
setwd("../..")
save(gpl.4291.M.avg, file="gpl.4291.M.avg.RData")
save(gpl.4293.M.avg, file="gpl.4293.M.avg.RData")
save(gpl.5774.M.avg, file="gpl.5774.M.avg.RData")

save(gpl.4291.rv.M, file="gpl.4291.rv.M.RData")
save(gpl.4293.rv.M, file="gpl.4293.rv.M.RData")
save(gpl.5774.rv.M, file="gpl.5774.rv.M.RData")
