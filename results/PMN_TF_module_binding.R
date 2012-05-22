# Script to compute to how many module members each assigned TF binds.
# Also, calculate the distribution of edge scores
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

#Functions
get_module <- function(modID, modules){
  return(modules[modules$moduleID==modID,2])
}

get_tf_weight <- function(tf, pDNA){
  tf.weight <- pDNA[pDNA$tf==tf,]$weight[1]
  return(tf.weight)
}

get_tf <- function(modID, modules, pDNA){
  mod.members <- get_module(modID, modules)
  mod.tfs <- table(pDNA[pDNA$target%in%mod.members,]$tf)
  mod.tfs.df <- as.data.frame(mod.tfs)
  if (nrow(mod.tfs.df)==0){
    return(data.frame(Var1=c(NA), Freq=c(NA), modID=c(modID)))
  }
  mod.tfs.df$moduleID <- modID
  return(mod.tfs.df)
}

get_module_tfs <- function(modIDs, modules, pDNA){
  results <- get_tf(modIDs[1], modules, pDNA)
  for (modID in modIDs[2:length(modIDs)]){
    results <- rbind(results, get_tf(modID, modules, pDNA))
  }
  
  return(results)
}

parse_pathways <- function(path.raw){
  pathways <- list()
  
  for (i in c(1:nrow(path.raw))){
    modID <- path.raw[i,1]
    pathway <- unlist(strsplit(path.raw[i,2], "\t", fixed=T))
    #Each parent has its own pathway, append the parent to the ID
    modID <- paste(c(modID, "_", pathway[1]), collapse="")
    pathways[[modID]] <- pathway
    
  }
  
  return(pathways)
}

tf_binding <- function(tf, module, pDNA){
  tf.targets <- pDNA[pDNA$tf==tf,2]
  mod.bound <- sum(module%in%tf.targets)
  return(mod.bound)
}

calc_tf <- function(pathways, modules, pDNA){
  
  tf.stats <- data.frame(modID=c(), parent=c(), tf=c(), 
                         n_bound=c(), perc_bound=c())
  for (pathID in names(pathways)){
    parsed <- unlist(strsplit(pathID, "_", fixed=T))
    modID <- parsed[1]
    parent <- parsed[2]
    
    pathway <- pathways[[pathID]]
    tf <- pathway[length(pathway)]
    
    module <- get_module(modID, modules)
    n_bound <- tf_binding(tf, module, pDNA)
    perc_bound <- n_bound / length(module)
    tf.weight <- get_tf_weight(tf, pDNA)
    
    tf.stats <- rbind(tf.stats, c(modID, parent, tf, tf.weight, 
                                  n_bound, perc_bound))
  }
  colnames(tf.stats) <- c("moduleID", "parent", "tf", "tf.weight",
                          "n_bound", "perc_bound")
  return(tf.stats)
}

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

pathways <- parse_pathways(path.raw)



tf.stats <- calc_tf(pathways, modules, tfs)
module.tfs <- get_module_tfs(good.modules, modules, tfs)

#mod2 TFs
mod2.tfs <- get_tf("mod2", modules, tfs)
write.table(mod2.tfs, "data/results/PMN_mod2_TFs.txt",
            col.names=T, row.names=F, quote=F, sep="\t")

write.table(tf.stats, "data/results/PMN_TF_module_binding.txt",
            row.names=F, col.names=T, quote=F, sep="\t")

