#Script to map the RvIDs to common names for the module membership and parents and pathways

#Written by Joe Reistetter

options(stringsAsFactor=F)
library(biomaRt)

get_parents <- function(df){
  modIDs <- unique(df$moduleID)
  
  out <- data.frame(moduleID=vector(mode="character"),
                    parent=vector(mode="character"))
  
  for (mod in modIDs){
    mod.rows <- df[df$moduleID==mod,]
    mod.parents <- mod.rows[1,2]
    
    #Check how many parents
    if (nchar(mod.parents) < 6){
      out <- rbind(out, c(mod, mod.parents))
    }
    else{
      mod.parents <- unlist(strsplit(mod.parents, " ", fixed=T))
      for (par in mod.parents){
        out <- rbind(out, c(mod, par))
      }
    }
  }
  colnames(out) <- c("moduleID", "gene")
  return(out)
}
#Load data
setwd("~/Dropbox/thesis_work/")
    
modules.raw <- read.table("PMN_output/4.17.30_mods_members.txt",
      head=T, sep="\t")

good.modules <- read.table("data/results/PMN_good_modules.txt",
                           head=T, sep="\t")

module_members <- modules.raw[modules.raw$moduleID%in% good.modules$moduleID,]
module_members$parent <- FALSE

#Extract parents and make a new DF with parents included.

parents.raw <- read.table("PMN_output/4.17.30_mods_parsed.txt",
                          head=T, sep="\t", stringsAsFactors=F)

parents.raw <- parents.raw[parents.raw$moduleID %in% good.modules$moduleID,]
parents <- get_parents(parents.raw)
parents$parent <- TRUE

modules <- rbind(module_members, parents)

#annotate data frame to gene names


bac <- useMart('bacteria_mart_13', dataset='myc_30_gene')

gene.ids.1 <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                    filters="external_gene_id", 
                    values=unique(gene.names.1), mart=bac)