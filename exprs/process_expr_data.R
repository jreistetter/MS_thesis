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
setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/exprs")
#My laptop wd
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/GSE8786/")


##########################

#  GSE8786

##########################

#First, need to create Targets frame for use in the read.maimages func

#Pull info from the GEO repository so we can get annotations for each array
dir.create("./GS8786")
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

##Now that the data is read in, do some QA/QC
dir.create("../QA/prenormMA")
plotMA3by2(g.8786.rg, path="../QA/prenormMA")

#Normalize the arrays
g.8786.bc <- backgroundCorrect(g.8786.rg, method="normexp", offset=50)
g.8786.MA <- normalizeWithinArrays(g.8786.rg, method="loess")
g.8786.MAq <- normalizeBetweenArrays(g.8786.rg, method="Aquantile")


#QA plots to see if normalization worked
setwd("./QA")
png("uncorrected_densities.png")
plotDensities(g.8786.rg)
dev.off()

png("corrected_densities.png")
plotDensities(g.8786.MA)
dev.off()

plotDensities(g.8786.MAq)

dir.create("./postnormMA")
plotMA3by2(g.8786.MA, path="./postnormMA")



#GEO IDs for the five datasets
GEO.ids <- c("GSE8839", "GSE16146", "GSE10391", "GSE8786", "GSE9331")

#Download the supplementary data (raw files) for each dataset
i=1
dl.file.info <- list()
for(i in 1:length(GEO.ids)){
  dl.file.info[[GEO.ids[i]]] <- getGEOSuppFiles(GEO.ids[i], makeDirectory=F)
}

#This is returning a workable file, need to use the GEO tutorial to extract it
g.16146 <- getGEO(GEO.ids[2], destdir="./", GSEMatrix=F)
Table(GSMList(g.16146)[[1]])[1:5,]

#Probably will loop through each element of GSMList and write out to text file,
#then read back in (or maybe just make into new Dfs)

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
plat1.assay <- assayData(gse.plat1)








probesets <- as.character(Table(GSMList(gse.plat1)[[1]])$ID_REF)
