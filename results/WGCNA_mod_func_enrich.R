# Script to analyze the modules obtained from WGCNA analysis for overrepresentation
# of GO terms to assign biological function to modules.
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

library(GOstats)

mod.hyper <- function(name, mod_genes, universe_genes, 
                      ontology, pvalue=0.05, categorySize=5){
  params <- GSEAGOHyperGParams(name=name,
                               geneSetCollection=H37Rv.gsc,
                               geneIds=mod_genes,
                               universeGeneIds=universe_genes,
                               ontology=ontology,
                               pvalueCutoff=0.05,
                               conditional=F,
                               testDirection="over")
  
  mod.test <- hyperGTest(params)
  results <- summary(mod.test, pvalue=pvalue, categorySize=categorySize)
  return(results)
}

mod.hyper.batch <- function(mod_members, universe, ontology, 
                            pvalue=0.05, categorySize=5){
  
  mods <- unique(mod_members$moduleID)
  mod1.genes <- mod_members[mod_members$moduleID==mods[1],2]
  results <- mod.hyper(mods[1], mod1.genes, universe,
                       ontology, pvalue, categorySize)
  
  #If no results returned, then put in NA values for the module
  if (nrow(results) == 0){
    results[1,] <- rep(NA,7)
  }
  results$moduleID <- mods[1]
  
  for (mod in mods[2:length(mods)]){
    mod.genes <- mod_members[mod_members$moduleID==mod,2]
    mod.results <- mod.hyper(mod, mod.genes, universe, 
                             ontology, pvalue, categorySize)
    
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

#Load the WGCNA modules
load("exprs/filt_pt5.net.RData")








