# Script to do network analysis of PMN PPI and pDNA network

# Actually need PPI network not pDNA.

options(stringsAsFactors=F)
library(igraph)

# Load data and process for igraph
setwd("~/Dropbox/thesis_work/data/")
pDNA <- read.table("protein-DNA/H37Rv.pdna.list", head=F, sep="\t")
PPI <- read.table("PPI/mtb.pp.list", head=F, sep="\t")

to_edge_list <- function(row){
  "Strips off the RV and optional letter at end for edge format"
  row <- as.character(row)
  v1 <- substr(row[1], 3, 6)
  v2 <- substr(row[2], 3, 6)
  return(paste(v1, v2, sep=" "))
}

get_shortest <- function(graph, v1, v2){
  "Returns the length of the shortest path between two proteins"
  
  len <- length(unlist(get.shortest.paths(graph, v1, v2)))
  return(len)
}

edge_list <- apply(pDNA, 1, to_edge_list)
ppi.edge_list <- apply(PPI, 1, to_edge_list)

write.table(edge_list, "protein-DNA/H37Rv.edge_list", 
            row.names=F, col.names=F, quote=F)
write.table(ppi.edge_list, "PPI/H37Rv.ppi.edge_list", 
            row.names=F, col.names=F, quote=F)

# Load data into igraph

pDNA.graph <- read.graph("protein-DNA/H37Rv.edge_list.test", 
                         format="edgelist", directed=F)

PPI.graph <- read.graph("PPI/H37Rv.ppi.edge_list",
                        format="edgelist", directed=F)


# Explore module 13
# Parents: 3103 and 0792
# Assigned TF: 3246C
# 2711 and 3676 both have more targets in mod

get_shortest(PPI.graph, 3103, 3246)
#[1] 4
get_shortest(PPI.graph, 3103, 2711)
# [1] 5
get_shortest(PPI.graph, 3103, 3676)
#[1] 4

get_shortest(PPI.graph, 0792, 3246)
# [1] 0
# Warning message:
#   In get.shortest.paths(graph, v1, v2) :
#   At structural_properties.c:917 :Couldn't reach some vertices
get_shortest(PPI.graph, 0792, 2711)
# [1] 0
# Warning message:
#   In get.shortest.paths(graph, v1, v2) :
#   At structural_properties.c:917 :Couldn't reach some vertices
get_shortest(PPI.graph, 0792, 3676)
# [1] 0
# Warning message:
#   In get.shortest.paths(graph, v1, v2) :
#   At structural_properties.c:917 :Couldn't reach some vertices