# Script to compute to how many module members each assigned TF binds.
# Also, calculate the distribution of edge scores
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

setwd("~/Dropbox/thesis_work/")

#########################
#
#   Load TF binding data and modules
#
#########################
tfs <- read.table("PMN_data/4.12.2012/H37Rv.pdna.list",
                  head=F, sep="\t", 
                  col.names=c("tf", "target", "weight", "direction"))

#Load module membership and pathways
good.modules <- read.table("data/results/PMN_good_modules.txt",
                           head=T, sep="\t")[,1]

modules.raw <- read.table("PMN_output/4.17.30_mods_members.txt",
                          head=T, sep="\t")

modules <- modules.raw[modules.raw$moduleID %in% good.modules,]

write.table(modules, "PMN_output/good_modules.txt",
            row.names=F, col.names=T, quote=F, sep="\t")

#########################
#
#   Load pathways
#
#########################

path.raw <- read.table("PMN_output/4.17.30_mods_pathways_for_annotation.txt",
                       head=F, sep=";")

path.raw <- path.raw[path.raw[,1] %in% good.modules,]

pathways <- list()

for (i in c(1:nrow(path.raw))){
  modID <- path.raw[i,1]
  pathway <- unlist(strsplit(path.raw[i,2], "\t", fixed=T))
  #Each parent has its own pathway, append the parent to the ID
  modID <- paste(c(modID, "_", pathway[1]), collapse="")
  pathways[[modID]] <- pathway
  
}


tf_binding <- function(tf, module, pDNA){
  tf.targets <- pDNA[pDNA$tf==tf,2]
  mod.bound <- sum(module%in%tf.targets)
  return(mod.bound)
}



