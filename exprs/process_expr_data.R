#Script to download the 5 GEO datasets that will be used to build the PMN.

library(GEOquery)

setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs")

#GEO Numbers
geos <- c("GSE8839", "GSE16146", "GSE10391", "GSE8786", "GSE9331")

gse8839 <- getGEO(geos[1], destdir=".")