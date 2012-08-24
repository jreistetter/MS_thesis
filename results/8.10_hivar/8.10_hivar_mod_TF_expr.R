# Script to look at characteristics of PMN modules for 8.10_hivar including:
#   TF binding
#   Overall module expression

options(stringsAsFactors=F)

#Functions
parse_TF <- function(parse_df){
  #   Parses out the TF for each module from the pathways parsed
  #   by 8.10.12_30_mods_0.1_highvar_parse.py.
  #   
  #   Params:
  #     parse_df: df read in from the parsed txt file, split on ":". First column 
  #               is the module ID, the second column is a tab delimted char vector
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

#Load data
setwd("~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/")

modules <- read.table("output/8.10_hivar_modules.txt",
                      sep="\t", head=T)

pathways <- read.table("8.10.12_30_mods_0.1_highvar_mods_pathways.txt",
                       sep=":", head=F)