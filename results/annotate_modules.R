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


##############################
#
#   Annotate the modules
#
##############################
    
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

bac <- useMart('bacteria_mart_13', dataset='myc_30_gene')

gene.ids <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                    filters="tuberculist", 
                    values=modules$gene, mart=bac)

gene.ids$tuberculist <- toupper(gene.ids$tuberculist)
modules.annotated <- merge(modules, gene.ids, by.x="gene", by.y="tuberculist", all.x=T)[,c(2:4)]
colnames(modules.annotated) <- c("module", "parent", "gene")

modules.annotated <- modules.annotated[with(modules.annotated, order(module, -parent)),]

write.table(modules.annotated, "data/results/modules_annotated.txt",
            col.names=T, sep="\t", quote=F, row.names=F)


