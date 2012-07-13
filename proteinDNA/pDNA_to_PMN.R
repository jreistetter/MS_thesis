# Script to assign protein-DNA edges from regulator-target pairs
# and operon predictions.
# 
# Written by Joe Reistetter
#
# 
#   Results:
#
#     p-DNA edges between 869 genes
#
#     534 pairs from the literature
#     403 pairs from MycoRegNet
#     181 pairs from MtbRegList
#     _________
#     1,118 regulator/target pairs in total, 173 of which were duplicated
#     945 *unique* regulator/target pairs, final
#     
#
#     4560 total edges
#     3269 duplicated edges
#     1276 protein- DNAedges remain after removing duplicates and 15 auto-regulation edges
#     
# 
# Workflow:
#   Do this for each operon data set.
# 
#   For each regulator
#     For each target
#       1. look up target operon ID in operon_genes dataframe
#       2. get character vector of operon members
#       3. Get index of target in vector
#       4. If on + strand (not c):
#           Create edges between regulator and target
#           Create edges between regulator and vector idx > target idx
#          If on c strand:
#           Create edges between regulator and target
#           Create edges between regulator and vector idx < target idx
#   
#   After doing that for each operon data set, remove any duplicate edges.
#   At this time, too worried about non-independence of operon predictions to do
#   a confidence score based on the number of DBs edge is in.

#At OHSU Dropbox
root_path <-"~/Dropbox/thesis_work"


setwd(paste(root_path, "/data/protein-DNA/RData", sep=""))
source(paste(root_path, "/code/proteinDNA/pdna_funcs.R", sep=""))

#Import regulator/target pairs from databases
load("mtbreglist.pairs.RData")
nrow(mtbreglist.pairs)
#[1] 181 pairs from MtbRegList

load("myco.pairs.TFBS.RData")
load("myco.pairs.ortho.RData")
nrow(myco.pairs.TFBS) + nrow(myco.pairs.ortho)
#[1] 403 from MycoRegNet

#Import regulator/target pairs from literature
setwd("../literature/")
lit.reg_targ <- read.table("pDNA_literature.txt", head=T,
                           sep='\t', stringsAsFactors=F)

dim(lit.reg_targ)
#[1] 534   3
#534 pairs
#3rd column is the literature source, don't neeed that.

reg_target.raw <- rbind(mtbreglist.pairs, myco.pairs.TFBS)
reg_target.raw <- rbind(reg_target.raw, myco.pairs.ortho)
reg_target.raw <- rbind(reg_target.raw, lit.reg_targ[,c(1,2)])

nrow(reg_target.raw)
#[1] 1118 pairs

reg_target.raw <- as.data.frame(lapply(reg_target.raw, toupper),
                                stringsAsFactors=F)

#Remove duplicate edges, if any
reg_target.raw$pair <- paste(reg_target.raw$regulator, 
                         reg_target.raw$target, sep="")

sum(duplicated(reg_target.raw$pair))
#[1] 173

reg_target <- reg_target.raw[!duplicated(reg_target.raw$pair),]
nrow(reg_target)
#[1] 945 unique pairs

nrow(reg_target.raw) - nrow(reg_target) == 173 #True, got all dupes

###############################
#
#      DOOR
#
###############################
setwd("../RData/")


load("door.op.list.RData")
load("door_genes.RData")

#Convert to uppercase
door_genes <- as.data.frame(lapply(door_genes, toupper),
                            stringsAsFactors=F)

#Loop through regulator/targets and assign edges
door.edges <- assign_pDNA_edge(reg_target, door_genes, door.op.list)


###############################
#
#      Microbes Online
#
###############################

load("microbes.op.list.RData")
load("microbes_genes.RData")

microbes_genes <- as.data.frame(lapply(microbes_genes, toupper),
                            stringsAsFactors=F)

microbes_online.edges <- assign_pDNA_edge(reg_target, microbes_genes, microbes.op.list)


###############################
#
#      ODB
#
###############################

load("ODB.op.list.RData")
load("ODB_genes.RData")



ODB_genes <- as.data.frame(lapply(ODB_genes, toupper),
                           stringsAsFactors=F)

ODB.edges <- assign_pDNA_edge(reg_target, ODB_genes, ODB.op.list)


#Combine edges

pDNA.edges.raw <- rbind(door.edges, microbes_online.edges)
pDNA.edges.raw <- rbind(pDNA.edges.raw, ODB.edges)

pDNA.edges.df <- as.data.frame(pDNA.edges.raw, stringsAsFactors=F)

pDNA.edges.df$pair <- paste(pDNA.edges.df$protein,
                             pDNA.edges.df$gene, sep="")

pDNA.edges.df <- as.data.frame(lapply(pDNA.edges.df, toupper),
                               stringsAsFactors=F)

nrow(pDNA.edges.df)
#[1] 4560 edges

sum(duplicated(pDNA.edges.df$pair))
#[1] 3269


#Remove dupes
pDNA.edges <- pDNA.edges.df[!duplicated(pDNA.edges.df[,3]),]
nrow(pDNA.edges)
#[1] 1291

#Check that dupes removed:
nrow(pDNA.edges.raw) - nrow(pDNA.edges) == 3269
#[1] TRUE

#Remove any auto-regulation pairs
pDNA.edges <- pDNA.edges[-which(pDNA.edges[,1]==pDNA.edges[,2]),]
nrow(pDNA.edges)
#[1] 1276, so 15 were autoregulation

#remove pairs column
pDNA.edges <- pDNA.edges[,-3]


#Calculate confidence score based on PMN paper

edge.genes <- unique(c(pDNA.edges[,1], pDNA.edges[,2]))
length(edge.genes)
#[1] 869 genes

#initialize list to hold degrees
node_degrees <- vector("list", length(edge.genes))
names(node_degrees) <- edge.genes
node_degrees <- lapply(node_degrees, function(x) x <- 0)

#Add 1 to the degree of each protein in an edge for STRING
for (node in c(pDNA.edges[,1], pDNA.edges[,2])){
  node_degrees[[node]] <- as.integer(node_degrees[[node]]) + 1
}


#From paper, the confidence score for an edge is given as:
#
#  p = 1 / 2 + exp(-0.2x + 5)
#
#

calc_conf <- function(node_degree){
  p <- 1 / (2 + exp(-0.2*node_degree + 5))
  return(p)
}

pDNA.edges$conf <- 0
add_degrees <- function(df, node_degrees){
  for (i in c(1:nrow(df))){
    degree <- node_degrees[[df[i,1]]]
    df[i,3] <- calc_conf(degree)
  }
  
  return(df)
}

pDNA.edges <- add_degrees(pDNA.edges, node_degrees)

#Assign forward direction to edge, see PMN docs
pDNA.edges$direction <- 1

save(pDNA.edges, file="pDNA.edges.RData")

write.table(pDNA.edges, file="../H37Rv.pdna.list", quote=F, sep="\t", row.names=F,
            col.names=F)


#Write out the transcription factors file, which are all regulators in the pDNA network.
tfs <- unique(pDNA.edges$protein)

write.table(tfs, file="../H37Rv.tfs", quote=F, sep="\t", row.names=F, col.names=F)



