# Script to analyze the modules obtained from 4.17 run for overrepresentation
# of GO terms to assign biological function to modules.
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

mod.hyper <- function(name, mod_genes, universe_genes, ontology){
  params <- GSEAGOHyperGParams(name=name,
                               geneSetCollection=H37Rv.gsc,
                               geneIds=mod_genes,
                               universeGeneIds=universe_genes,
                               ontology=ontology,
                               pvalueCutoff=0.05,
                               conditional=F,
                               testDirection="over")
  
  mod.test <- hyperGTest(params)
  return(mod.test)
}

#########################
#     Main
#########################

#Load data

setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data")

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

my.universe <- mod_members$gene


#Example module
mod23 <- mod_members[mod_members$moduleID=="mod23",2]
mod23.BP <- mod.hyper("mod23BF", mod23, universe, "BP")
mod23.BP.2 <- mod.hyper("mod23BF2", mod23, my.universe, "BP")
head(summary(mod23.BP.2))

mod23.CC <- mod.hyper("mod23BF", mod23, universe, "CC")
head(summary(mod23.CC))

