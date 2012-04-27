# Script to import gene-GO term mappings obtained from AmiGO and create a
# database object for use in topGO.
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

#Functions

get_rv <- function(pipe_list){
  #Extracts Rv ID from pipe-delimted list
  ids <- unlist(strsplit(pipe_list, '|', fixed=T))
  rv_id <- ids[grepl("Rv", ids, fixed=T)]
  if (length(rv_id) == 0){
    return(NA)
  }
  if (length(rv_id) > 1){
    return(rv_id[2])
  }
  return(rv_id)
}

#########################
#       Main
#########################
library(topGO)
#Load data

#OHSU
setwd("~/Dropbox/thesis_work/data/GO")

go.raw <- read.table("amigo_all_genes.txt", quote="",
                     head=F, sep='\t', comment.char="",)

dim(go.raw)
#[1] 7343   15

#The 5th column contains the GO term. 11th column contains pipe
#delimted text of different gene IDs, one of which is an Rv ID that we want.
#Not all rows have genes that are in H37Rv, so remove those.

go.rv <- go.raw[grepl("Rv", go.raw[,11], fixed=T),]

dim(go.rv)
#[1] 7224   15

#make new object of just Rv ID and GO mappings
#Use get_rv to extract the Rv ID from the pipe list

rv.ids <- toupper(sapply(go.rv[,11], get_rv, USE.NAMES=F))

#Check filtering
stopifnot(length(rv.ids) == nrow(go.rv)) #Make sure that no rows removed
stopifnot(sum(grepl("RV", rv.ids, fixed=T))==nrow(go.rv)) #Check that all are Rvs
stopifnot(sum(is.na(rv.ids))==0) #check for NAs

#Create data frame with GO term and gene
go.mtb.mappings <- data.frame(GO=go.rv[,5], gene=rv.ids)

stopifnot(sum(is.na(go.mtb.mappings$GO))==0)


# Write out mappings in format to read back in with topGO
#Format: geneID\tgo1,go2,go3....

#Build list of genes and their associated GO terms
geneGO <- list()
gene.ids <- unique(go.mtb.mappings$gene)
length(gene.ids)
#[1] 2312

out_file <- file("GO_mappings.txt", "w")

for (gene in gene.ids){
  gene.gos <- go.mtb.mappings[go.mtb.mappings$gene == gene,1]
  gos.csv <- paste(gene.gos, collapse=", ")
  line <- paste(c(gene, "\t ", gos.csv), collapse="")
  write(line, out_file, append=T)
}

close(out_file)

#Convert to GO object
mtb.GO.db <- readMappings("GO_mappings.txt")

#Save for use in later script
save(mtb.GO.db, file="mtb.GO.db.RData")




