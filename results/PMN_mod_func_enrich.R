# Script to analyze the modules obtained from 4.17 run for overrepresentation
# of GO terms to assign biological function to modules.
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

library(GOstats)

mod.hyper <- function(name, mod_genes, universe_genes, 
                      ontology, pvalue=0.05, categorySize=5, conditional=F){
  params <- GSEAGOHyperGParams(name=name,
                               geneSetCollection=H37Rv.gsc,
                               geneIds=mod_genes,
                               universeGeneIds=universe_genes,
                               ontology=ontology,
                               pvalueCutoff=0.05,
                               conditional=conditional,
                               testDirection="over")
  
  mod.test <- hyperGTest(params)
  results <- summary(mod.test, pvalue=pvalue, categorySize=categorySize)
  return(results)
}

mod.hyper.batch <- function(mod_members, mod_parents, universe, ontology, 
                            pvalue=0.05, categorySize=5, conditional=F){
  
  mods <- unique(mod_members$moduleID)
  mod1.genes <- mod_members[mod_members$moduleID==mods[1],2]
  mod1.parents <- mod_parents[mod_parents$moduleID==mods[1],2]
  mod1.parents <- unique(unlist(strsplit(mod1.parents, " ", fixed=T)))
  mod1.genes <- c(mod1.genes, mod1.parents)
  
  results <- mod.hyper(mods[1], mod1.genes, universe,
                       ontology, pvalue, categorySize, conditional)
  
  #If no results returned, then put in NA values for the module
  if (nrow(results) == 0){
    results[1,] <- rep(NA,7)
  }
  results$moduleID <- mods[1]
  
  for (mod in mods[2:length(mods)]){
    mod.genes <- mod_members[mod_members$moduleID==mod,2]
    mod.parents <- mod_parents[mod_parents$moduleID==mod,2]
    mod.parents <- unique(unlist(strsplit(mod.parents, " ", fixed=T)))
    mod.genes <- c(mod.genes, mod.parents)
    mod.results <- mod.hyper(mod, mod.genes, universe, 
                             ontology, pvalue, categorySize, conditional)
    
    if (nrow(mod.results) == 0){
      mod.results[1,] <- rep(NA,7)
    }
    
    mod.results$moduleID <- mod
    results <- rbind(results, mod.results)
  }
  results$ontology <- ontology
  
  return(results[complete.cases(results),])
}

#########################
#     Main
#########################

#Load data

#OHSU
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data")

#Laptop
setwd("~/schoolDB/Dropbox/thesis_work/data")

#Load the GO object and universe of genes with GO mappings
load("GO/H37Rv.gsc.RData")
load("GO/universe.RData")

#Module membership
mod_members.raw <- read.table("../PMN_output/4.17.30_mods_members.txt",
                              head=T, sep='\t')
dim(mod_members.raw)
#[1] 3458    2

#Filter out members not in the universe:
mod_members <- mod_members.raw[mod_members.raw$gene %in% universe,]
dim(mod_members)
#[1] 2075    2, so ~1400 have no associated GO term

mod_parents <- read.table("../PMN_output/4.17.30_mods_parsed.txt",
                          head=T, sep='\t')

#Set the universe to be all genes present in the analysis that have GO terms
my.universe <- mod_members$gene

#Only analyze modules with at least 1 meaningful probability

#Load in module statistics
mod.stats <- read.table("../PMN_output/4.17_30mods_genes_pathsizes.txt",
                        head=T, sep='\t')
#Filter on number of probabilities > 0.4
mod.good <- mod.stats[mod.stats$thresh.0.2 > 0,]$moduleID
mod_members.good <- mod_members[mod_members$moduleID %in% mod.good,]

#Not conditioned
all.mods.BP <- mod.hyper.batch(mod_members.good, mod_parents, my.universe, "BP")
all.mods.CC <- mod.hyper.batch(mod_members.good, mod_parents, my.universe, "CC")
all.mods.MF <- mod.hyper.batch(mod_members.good, mod_parents, my.universe, "MF")

colnames(all.mods.BP)[1] <- "GO_ID"
colnames(all.mods.CC)[1] <- "GO_ID"
colnames(all.mods.MF)[1] <- "GO_ID"

mods.GO.enrichment <- rbind(all.mods.BP, all.mods.CC)
mods.GO.enrichment <- rbind(mods.GO.enrichment, all.mods.MF)

write.table(mods.GO.enrichment, file="./GO/4.17_module_GO_enrichment.txt",
            col.names=T, row.names=F, quote=F, sep='\t')



#Conditioned on GO structure
cond.mods.BP <- mod.hyper.batch(mod_members.good, mod_parents,
                                my.universe, "BP", conditional=T)
cond.mods.CC <- mod.hyper.batch(mod_members.good, mod_parents,
                                my.universe, "CC", conditional=T)
cond.mods.MF <- mod.hyper.batch(mod_members.good, mod_parents,
                                my.universe, "MF", conditional=T)

colnames(cond.mods.BP)[1] <- "GO_ID"
colnames(cond.mods.CC)[1] <- "GO_ID"
colnames(cond.mods.MF)[1] <- "GO_ID"

mods.cond.GO.enrichment <- rbind(cond.mods.BP, cond.mods.CC)
mods.cond.GO.enrichment <- rbind(mods.cond.GO.enrichment, cond.mods.MF)

write.table(mods.cond.GO.enrichment, file="GO/4.17module_func_enrichment_conditional.txt",
            col.names=T, row.names=F, quote=F, sep='\t')
