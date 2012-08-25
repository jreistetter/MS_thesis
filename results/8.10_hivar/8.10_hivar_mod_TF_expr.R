# Script to look at characteristics of PMN modules for 8.10_hivar including:
#   TF binding assigned by PMN
#   pDNA interactions in each module
#   Overall module expression

options(stringsAsFactors=F)

#Functions

setwd("~/Dropbox/thesis_work/code/")
source("thesis_funcs.R")
parse_TF <- function(parse_df){
  #   Parses out the TF for each module from the pathways parsed
  #   by 8.10.12_30_mods_0.1_highvar_parse.py.
  #   
  #   Params:
  #     parse_df: df read in from the parsed txt file, split on ":". First column 
  #             is the module ID, the second column is a tab delimted char vector
  #               of the protein pathway from module parent to TF (TF is last).
  
  out_df <- data.frame(modID=c(), tf=c())
  for (i in 1:nrow(parse_df)){
    modID <- parse_df[i,1]
    path.raw <- parse_df[i,2]
    path_members <- unlist(strsplit(path.raw, "\t"))
    tf <- path_members[length(path_members)]
    out_df <- rbind(out_df, c(modID, tf))
    }
  colnames(out_df) <- c("modID", "tf")
  return(out_df)
}

write_mod_binding <- function(modID, mod.df, pDNA){
#   Retrieves the known TF binding for a module.
#   Writes out assigned TF and TF binding table to file
#     
#     Params:
#       modID - module ID (str)
#       mod.df - modules df
#       pDNA - df of pDNA edges
    
  mod.TFs <- get_tf(modID, modules, pDNA)
  mod.assigned.TF <- unique(modules_tf[modules_tf$modID==modID,2])
  fname <- paste(modID, "_tf_binding.txt", sep="")
  out_f <- file(fname, "w")
  write("assigned TFs:", out_f)
  write(mod.assigned.TF, out_f, sep="\t")
  write("\n", out_f)
  write.table(mod.TFs, out_f, append=T, quote=F, sep="\t",
              row.names=F, col.names=F)
  close(out_f)
}

####################
#       Main
#
####################

#Load data
setwd("~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/")

modules <- read.table("output/8.10_hivar_modules.txt",
                      sep="\t", head=T)

pathways <- read.table("8.10.12_30_mods_0.1_highvar_mods_pathways.txt",
                       sep=":", head=F)

setwd("~/Dropbox/thesis_work/data/")

pDNA.raw <- read.table("protein-DNA/H37Rv.pdna.list", sep="\t")
pDNA <- pDNA.raw[,1:2]
colnames(pDNA) <- c("tf", "target")

modules_tf <- parse_TF(pathways)

write.table(modules_tf, 
            "~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/output/8.10_hivar_module_tfs.txt",
            quote=F, sep="\t", row.names=F, col.names=T)

modIDs <- unique(modules$moduleID)


#For each module, write out table of TFs targetting module and also the assigned
# TF

setwd("~/Dropbox/thesis_work/data/8.10_highvar_results/Module_TF_binding/")

for (ID in modIDs){
  write_mod_binding(ID, modules, pDNA)
}





