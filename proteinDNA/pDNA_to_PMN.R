# Script to assign protein-DNA edges from regulator-target pairs
# and operon predictions.
# 
# Written by Joe Reistetter
#
# 
#   Results:
#     434 regulator/target pairs in total, 150 of which were duplicated
#     728 protein-DNA edges written out to H37Rv.pdna.list
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
root_path <-"/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work"
#My laptop Dropbox
#root_path <- "~/schoolDB/Dropbox/thesis_work"

setwd(paste(root_path, "/data/protein-DNA/RData", sep=""))
source(paste(root_path, "/code/proteinDNA/pdna_funcs.R", sep=""))

#Import regulator/target pairs
load("mtbreglist.pairs.RData")
load("myco.pairs.TFBS.RData")
load("myco.pairs.ortho.RData")

reg_target.raw <- rbind(mtbreglist.pairs, myco.pairs.TFBS)
reg_target.raw <- rbind(reg_target.raw, myco.pairs.ortho)

nrow(reg_target.raw) #584 pairs

reg_target.raw <- as.data.frame(lapply(reg_target.raw, toupper),
                                stringsAsFactors=F)

#Remove duplicate edges, if any
reg_target.raw$pair <- paste(reg_target.raw$regulator, 
                         reg_target.raw$target, sep="")

sum(duplicated(reg_target.raw$pair))
#150

reg_target <- reg_target.raw[!duplicated(reg_target.raw$pair),]

nrow(reg_target.raw) - nrow(reg_target) == 150 #True, got all dupes

###############################
#
#      DOOR
#
###############################



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

sum(duplicated(pDNA.edges.raw[,2]))
#[1] 1313, most likely reg/targs without an operon


#Remove dupes
pDNA.edges <- pDNA.edges.raw[!duplicated(pDNA.edges.raw[,2]),]

nrow(pDNA.edges.raw) - nrow(pDNA.edges) == 1313

save(pDNA.edges, file="pDNA.edges.RData")

write.table(pDNA.edges, file="../H37Rv.pdna.list", quote=F, sep="\t", row.names=F,
            col.names=F)

