# Script to assign protein-DNA edges from regulator-target pairs
# and operon predictions.
# 
# Written by Joe Reistetter
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
#root_path <-"/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work"
#My laptop Dropbox
root_path <- "~/schoolDB/Dropbox/thesis_work"

setwd(paste(root_path, "/data/protein-DNA/RData", sep=""))
source(paste(root_path, "/code/proteinDNA/pdna_funcs.R", sep=""))

#Import regulator/target pairs
load("mtbreglist.pairs.RData")
load("myco.pairs.TFBS.RData")
load("myco.pairs.ortho.RData")

reg_target <- rbind(mtbreglist.pairs, myco.pairs.TFBS)
reg_target <- rbind(reg_target, myco.pairs.ortho)

nrow(reg_target) #584 pairs




