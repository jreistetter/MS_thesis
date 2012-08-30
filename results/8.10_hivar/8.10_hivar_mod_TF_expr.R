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

write_mod_binding <- function(modID, mod.df, pDNA, modules.tf){
#   Retrieves the known TF binding for a module.
#   Writes out assigned TF and TF binding table to file
#     
#     Params:
#       modID - module ID (str)
#       mod.df - modules df
#       pDNA - df of pDNA edges
#       modules.tf - df of TFs assigned to modules, from parse_TF()
  
  module <- get_module(modID, mod.df)
  mod.TFs <- get_tf(modID, modules, pDNA)
  mod.assigned.TF <- unique(modules_tf[modules_tf$modID==modID,2])
  fname <- paste(modID, "_tf_binding.txt", sep="")
  out_f <- file(fname, "w")
  write(paste("assigned TF:\t", mod.assigned.TF), out_f)
  write(paste("size:\t", length(module)), out_f)
  write.table(mod.TFs[,1:2], out_f, append=T, quote=F, sep="\t",
              row.names=F, col.names=F)
  
  n_no_mod <- sum(!(module %in% pDNA$target))
  write(paste("None assigned", n_no_mod, sep="\t"), out_f)
  close(out_f)
}

check_binding <- function(modID, modules, module_tfs, pDNA){
  #   Check if PMN assigned TF has most targets in assigned module
  #   Params:
  #     modID - module ID to be checked (str)
  #     modules - modules df
  #     module_tfs - df of assigned module TFs returned by parse_TF()
  #     pDNA - df of pDNA edges
  #     
  #     Returns:
  #       TRUE - assigned tf has most targets in module
  #       FALSE - assigned tf does not have most targets in module
  
  mod.tfs <- get_tf(modID, modules, pDNA)
  assigned <- unique(module_tfs[module_tfs$modID==modID,])
  
  max.tf.idx <- which.max(mod.tfs$Freq)
  max.tf <- mod.tfs[max.tf.idx,1]
  return(max.tf %in% assigned)
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
# Next step is for each module check that the assigned TF in PPI is
# the TF with the most targets in the module. i.e. the max of the Freq col
# in the df returned by get_tf(). Either keep a boolean value in a df for each
# module or something.

#Load in the stat sheet for modules and add column for TF assignment
mod.stats <- read.table(
  "~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/output/8.10.12_30_mods_0.1_highvar_genes_pathsizes.txt",
  head=T, sep="\t")

stat.ids <- mod.stats$moduleID
mod.stats$correct_tf_assigned <- sapply(mod.stats$moduleID,
                      function(id, modules, modules_tf, pDNA){
                      check_binding(id, modules, modules_tf, pDNA)}, 
                      modules=modules, modules_tf=modules_tf, pDNA=pDNA
)

mod.stats$correct_tf_assigned <- as.character(mod.stats$correct_tf_assigned)

mod.tf.assignments <- modules_tf[!duplicated(modules_tf$moduleID),]
colnames(mod.tf.assignments)[1] <- "moduleID"
mod.stats <- merge(mod.stats, mod.tf.assignments)

write.table(mod.stats,
            file="~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/output/8.10.12_30_mods_0.1_highvar_genes_pathsizes.txt",
            row.names=F, col.names=T, quote=F, sep="\t")


# Check which modules the dosR regulon was assigned to

# Load regulon membership
setwd("~/Dropbox/thesis_work/data/")
dosR_regulon <- toupper(read.table("known_regulons/DosR.txt")[,1])
dosR_assignments <- modules[modules$gene %in% dosR_regulon,]
dosR_assignments <- dosR_assignments[with(dosR_assignments, order(moduleID)),]
write.table(dosR_assignments, "8.10_highvar_results/dosR_regulon_assignments.txt",
            row.names=F, col.names=T, quote=F, sep="\t")

dosR_counts <- as.data.frame(table(dosR_assignments$moduleID))
colnames(dosR_counts) <- c("moduleID", "n_dosR_regulon")
write.table(dosR_counts, "8.10_highvar_results/dosR_regulon_counts.txt",
            row.names=F, col.names=T, quote=F, sep="\t")


