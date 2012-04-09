#Script for GSE8839 that:
# 1 - download and import the raw intensities from 2-dye arrays
# 2 - QA on the raw intensities
# 3 - Background correct and normalize the intensities
# 4 - QA on corrected and normalized intensities
# 5 - Discretize ratios to {-1, 0 1} for the treatment arrays

# Written by Joe Reistetter

library(GEOquery)
library(limma)

#At OHSU Dropbox
setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
##My laptop Dropbox
#setwd("~/schoolDB/Dropbox/thesis_work")
source("./code/exprs/exprs_funcs.R")

#setwd("./data/exprs/GSE8839")

#Pull info from the GEO repository so we can get annotations for each platform

#Downloads the soft files
g.8839.geo <- getGEO("GSE8839", destdir="./data/exprs/GSE8839", GSEMatrix=F) 
g.8839.gsm <- GSMList(g.8839.geo)

g.8839.geo.matrix <- getGEO("GSE8839", destdir="./data/exprs/GSE8839")

#Function to pull out the data table and save as text
save_table <- function(gsm, out_dir){
  dat <- Table(gsm@dataTable)
  accession <- gsm@header$geo_accession
  fname <- paste(accession, ".txt", sep="")
  write.table(dat, file=paste(out_dir,fname,sep=""), 
              sep="\t", quote=F, col.names=T, row.names=F)
  
}

dir.create("./data/exprs/GSE8839/GSE8839_RAW")
#setwd("./GSE16146_RAW")

#Save out to text
lapply(g.8839.gsm, save_table, out_dir = "./data/exprs/GSE8839/GSE8839_RAW/")

##########################

#  GSE8839 - GPL8561

##########################
#First, need to create Targets frame for use in the read.maimages func

#Gets phenotype data
g.8561 <- g.8839.geo.matrix[['GSE8839-GPL8561_series_matrix.txt.gz']]
g.8561.pheno <- phenoData(g.8561)
g.8561.pdata <- pData(g.8561.pheno) #Dataframe of all data associated with arrays
g.8561.annot <- featureData(g.8561)@data

g.8561.fnames <- rownames(g.8561.pdata)
g.8561.fnames <- sapply(g.8561.fnames, 
                          function(x) {paste(x, ".txt", sep="")},
                          USE.NAMES=F
                          )

#Get the GEO accessions for each array to use as an ID
g.8561.arraynames <- rownames(g.8561.pdata)

#Get the treatments for each array
g.8561.treatment <- as.character(g.8561.pdata$source_name_ch2)

#All the data is in the correct format, build it into a dataframe
#that fits the specs for input into the read.maimages function.
g.8561.targets <- data.frame(FileName=g.8561.fnames,
                               Cy3 = rep("control", length(g.8561.fnames)),
                               Cy5 = g.8561.treatment, stringsAsFactors=F)

rownames(g.8561.targets) <- g.8561.arraynames

##Set the column IDs for reading each dataset into the limma object
#column names: Ch1 Intensity (Mean), Ch1 Background (Median), Ch2 Background (Median), Ch2 Intensity (Mean) -- ch1 = Cy3 = green, ch2 = Cy5 = red
g.8561.cols <- list(R="CH2I_MEAN", G="CH1I_MEAN",
                      Rb="CH2B_MEDIAN", Gb="CH1B_MEDIAN")

#Define weight function to exclude spots that have a "FLAG" column value < 0 in the
#raw data
flagged <- function(x) as.numeric(x[,"FLAG"] > -1)

g.8561.rg <- read.maimages(g.8561.targets,
                             annotation=c("ID_REF"),
                             columns=g.8561.cols,
                             wt.fun=flagged,
                             path="./data/exprs/GSE8839/GSE8839_RAW")

g.8561.gal <- readGAL("./data/exprs/GSE8839/SMD_print_498.gal")

colnames(g.8561.gal)[5] <- "Reporter.Name"
#Generate a Layout object for the background correction and normalization
g.8561.rg$genes <- merge(g.8561.annot, g.8561.rg$genes, 
                         by.x = "ID", by.y="ID_REF", all.y=T)
g.8561.rg$printer <- getLayout(g.8561.gal)

##Now that the data is read in, do some QA/QC by looking at MA plots
dir.create("./data/exprs/GSE8839/GPL8561")
dir.create("./data/exprs/GSE8839/GPL8561/QA")
dir.create("./data/exprs/GSE8839/GPL8561/QA/prenormMA")
plotMA3by2(g.8561.rg, path="./data/exprs/GSE8839/GPL8561/QA/prenormMA", 
           main=paste(g.8561.targets$FileName, g.8561.targets$Cy5, paste=" - "))


#Normalize the arrays, may have to remove some if the artifacts remain
g.8561.bc <- backgroundCorrect(g.8561.rg, method="normexp", offset=50)
g.8561.bc.norm <- normalizeWithinArrays(g.8561.bc, method="loess")

#Need to filter so that only RV genes are present
g.8561.rv.idx <- grepl(pattern="RV", x=g.8561.rg$genes$ORF, fixed=T)
sum(g.8561.rv.idx) #4,686 features represented
g.8561.bc.norm.rv <- g.8561.bc.norm[g.8561.rv.idx,]
length(unique(as.character(g.8561.bc.norm.rv$genes$ORF))) #3924 genes

#Redo the MA plots and see if artifacts disappear
dir.create("./data/exprs/GSE8839/GPL8561/QA/postnormMA_RVfiltered")
plotMA3by2(g.8561.bc.norm.rv, 
           path="./data/exprs/GSE8839/GPL8561/QA/postnormMA_RVfiltered", 
           main=paste(g.8561.targets$FileName, g.8561.targets$Cy5, paste=" - "))

#QA plots to see if normalization worked
png("uncorrected_densities.png")
plotDensities(g.8561.rg)
dev.off()

png("./data/exprs/GSE8839/GPL8561/QA/corrected_densities.png")
plotDensities(g.8561.bc.norm.rv)
dev.off()

##Extract log-2 expression

#override remove_bad_spots function because gene name is
#in different named column for GSE16146
remove_bad_spots <- function(ma_list){
  probe.weights <- ma_list$weights
  probe.weights[probe.weights == 0] <- NA
  cleaned <- probe.weights * ma_list$M
  rownames(cleaned) <- as.character(ma_list$genes$ORF)
  return(as.data.frame(cleaned))
  
}

g.8561.rv.M <- remove_bad_spots(g.8561.bc.norm.rv)
g.8561.rv.M$gene <- rownames(g.8561.rv.M)
g.8561.rv.M.avg <- avg_probes(g.8561.rv.M, unique(g.8561.bc.norm.rv$genes$ORF))

dim(g.8561.rv.M.avg)
#[1] 3924  114
dim(g.8561.rv.M)
#[1] 4686  115

#Check, right number of arrays and rows.

save(g.8561.rv.M, file="./data/exprs/GSE8839/g.8561.rv.M.RData")
save(g.8561.rv.M.avg, file="./data/exprs/GSE8839/g.8561.rv.M.avg.RData")


