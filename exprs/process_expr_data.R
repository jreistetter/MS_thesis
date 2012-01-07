#Script to download the 5 GEO datasets that will be used to build the PMN.

library(GEOquery)

#At OHSU wd
setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs")
#My laptop wd
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs")


#GEO Numbers
geos <- c("GSE8839", "GSE16146", "GSE10391", "GSE8786", "GSE9331")

gse8786.raw <- getGEO(filename="GSE8786_RAW.tar")

gse8839 <- getGEO(geos[1], destdir=".")
gse8839.supp <- getGEOSuppFiles(geos[1])
gse8839.raw <- getGEO(filename="GSE8839_RAW.tar")

#This particular GSE has two platforms, getGEO returned a list containing both of them.
#Separate them out to two variables

#GPL2787 
gse.plat1 <- gse8839[[1]]
#GPL5714
gse.plat2 <- gse8839[[2]]

#First platform (GPL2787)
#See what data is stored in the ExpressionSet
plat1.pheno <- phenoData(gse.plat1) #contains sampleIDs, info about experiment blank
plat1.sampleIDs <- sampleNames(plat1.pheno) #get the sample names

#Take a look at the feature data
plat1.feats <- featureData(gse.plat1) #annotated dframe containing info about probes
plat1.probeIDs <- featureNames(plat1.feats) #probeIDs
plat1.pdata <- pData(plat1.feats) #info about probeID, including RV#









probesets <- as.character(Table(GSMList(gse.plat1)[[1]])$ID_REF)
