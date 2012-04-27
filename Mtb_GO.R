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








