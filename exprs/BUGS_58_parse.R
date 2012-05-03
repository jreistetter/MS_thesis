# Script to parse the BUGS-58 data into a usable R structure for
# input to Hotellings T2

# Written by Joe Reistetter

options(stringsAsFactor=F)
setwd("~/Dropbox/thesis_work/data/exprs/")

# Functions

get_RV <- function(row){
  first_split <- unlist(strsplit(as.character(row[1]), "-", fixed=T))
  id <- unlist(strsplit(first_split, " ", fixed=T))
  id <- paste(c("RV", id[2]), collapse="")
  return(id)
}

readBUGS <- function(samples, dir="./"){
  #Read in first array to see how many genes
  first_f <- paste(c(dir, samples[1,]$filename), collapse="/")
  first_arr <- read.table(first_f, header=T, sep="\t", quote="")
  ngenes <- nrow(first_arr)
  
  #Extract the Rv ID for rownames: MtH37Rv-0001 (1A1)
  gene.ids <- apply(first_arr, 1, get_RV)
  genes.df <- data.frame(geneID=gene.ids, stringsAsFactors=F)
  genes.df$first <- first_arr[,2]
  colnames(genes.df)[2] <- samples[1,]$filename
    
  for (i in c(2:nrow(samples))){
    f <- paste(c(dir, samples[i,]$filename), collapse="/")
    arr <- read.table(f, header=T, sep="\t", quote="")
    arr$geneID <- apply(arr, 1, get_RV)
    genes.df <- merge(genes.df, arr[,2:3],
                      all.x=T, all.y=T)
    
    colnames(genes.df)[i+1] <- samples[i,]$filename
  }
  rownames(genes.df) <- genes.df$geneID
  genes.df[,1] <- NULL
  return(genes.df)
  
}

#########################
#     Main
#########################

#Load sample annotations
samp.annot.raw <- read.table(
  "EBUGS58/E-BUGS-58_analysed/BUGS-58_sample_annotations.txt",
                         head=F, sep='\t', 
                         col.names=c("ID1", "ID2", "Pheno1"),
                         colClasses=rep("character", 3))


# Create file names
# Sample file format: TB73_012.txt

make_filename <- function(left, right){
  return(paste(c(left, right), collapse="_"))
}

samp.annot.raw$filename <- apply(samp.annot.raw, 1, 
                                   function(row){
                                     f <- paste(c(row[1], row[2]), 
                                                collapse="_")
                                     f <- paste(c(f, ".txt"),
                                                collapse="")
                                      return(f)
                                    }
                                   )

# Make separate columns for different sample characteristics

samp.annot.raw$celltype <- ""
samp.annot.raw$time <- ""
samp.annot.raw$donor <- ""

# The raw phenotype data looks like: DC 1h D1 (A1)

for (i in c(1:nrow(samp.annot.raw))){
  split <- unlist(strsplit(samp.annot.raw[i,3], " ", fixed=T))
  samp.annot.raw[i,]$celltype <- split[1]
  samp.annot.raw[i,]$time <- split[2]
  samp.annot.raw[i,]$donor <- split[3]
}

# Set the Aerobic sample time to NA
samp.annot.raw[samp.annot.raw$celltype=="Aerobic",6] <- NA

# Extract the analyzed data and save RData object
BUGS58.samples <- samp.annot.raw[,4:7]

save(BUGS58.samples, file="EBUGS58/BUGS58.samples.RData")

#########################
#
#   Read in the expression arrays
#
#########################

BUGS58.arrays <- readBUGS(BUGS58.samples, dir="EBUGS58/E-BUGS-58_analysed")
rownames(BUGS58.arrays) <- toupper(rownames(BUGS58.arrays))
save(BUGS58.arrays, file="EBUGS58/EBUGS58.arrays.RData")




