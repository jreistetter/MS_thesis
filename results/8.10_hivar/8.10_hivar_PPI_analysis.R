# Script to do network analysis of PMN PPI and pDNA network

options(stringsAsFactors=F)
library(igraph)

# Load data and process for igraph
setwd("~/Dropbox/thesis_work/data/protein-DNA/")
pDNA <- read.table("H37Rv.pdna.list", head=F, sep="\t")

to_edge_list <- function(row){
  row <- as.character(row)
  v1 <- substr(row[1], 3, 6)
  v2 <- substr(row[2], 3, 6)
  return(paste(v1, v2, sep=" "))
}

edge_list <- apply(pDNA, 1, to_edge_list)

write.table(edge_list, "H37Rv.edge_list", row.names=F, col.names=F, quote=F)
