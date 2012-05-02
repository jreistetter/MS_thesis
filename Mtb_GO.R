# Script to import gene-GO term mappings obtained from AmiGO and create a
# database object for use in topGO.
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

#Functions

get_rv <- function(pipe_list){
  #Extracts Rv ID from pipe-delimted list
  ids <- unlist(strsplit(pipe_list, '|', fixed=T))
  rv_id <- ids[grepl("Rv", ids, fixed=T)]
  if (length(rv_id) == 0){
    return(NA)
  }
  if (length(rv_id) > 1){
    return(rv_id[2])
  }
  return(rv_id)
}

mod_to_goeast <- function(modID, mod_members, mod_parents){
  mod.genes <- mod_members[mod_members$moduleID==modID,]$gene
  mod.parents <- mod_parents[mod_parents$moduleID==modID,2]
  mod.parents <- unique(unlist(strsplit(mod.parents, " ", fixed=T)))
  fname = paste(modID, ".txt", sep="")
  write.table(c(mod.genes, mod.parents), fname, 
              quote=F, row.names=F, col.names=F)
}

mods_GOEAST <- function(mod_members, mod_parents){
  modIDs <- unique(mod_members$moduleID)
  for (mod in modIDs){
    mod_to_goeast(mod, mod_members, mod_parents)
  }
}

#########################
#       Main
#########################
library(GOstats)
library(GSEABase)
#Load data

#OHSU
setwd("~/Dropbox/thesis_work/data/GO")

go.raw <- read.table("amigo_all_genes.txt", quote="",
                     head=F, sep='\t', comment.char="",)

dim(go.raw)
#[1] 7343   15

#The 5th column contains the GO term. 11th column contains pipe
#delimted text of different gene IDs, one of which is an Rv ID that we want.
#Not all rows have genes that are in H37Rv, so remove those.

go.rv <- go.raw[grepl("Rv", go.raw[,11], fixed=T),]

dim(go.rv)
#[1] 7224   15

#make new object of just Rv ID and GO mappings
#Use get_rv to extract the Rv ID from the pipe list

rv.ids <- toupper(sapply(go.rv[,11], get_rv, USE.NAMES=F))

#Check filtering
stopifnot(length(rv.ids) == nrow(go.rv)) #Make sure that no rows removed
stopifnot(sum(grepl("RV", rv.ids, fixed=T))==nrow(go.rv)) #Check that all are Rvs
stopifnot(sum(is.na(rv.ids))==0) #check for NAs

#Create data frame with GO term and gene
go.mtb.mappings <- data.frame(go_id=go.rv[,5], evidence=go.rv[,7], gene_id=rv.ids)
stopifnot(sum(is.na(go.mtb.mappings$go_id))==0)

universe <- unique(go.mtb.mappings$gene_id)
save(universe, file="universe.RData")

goFrame <- GOFrame(go.mtb.mappings, organism="Mtb H37Rv")
goAllFrame <- GOAllFrame(goFrame)

H37Rv.gsc <- GeneSetCollection(goAllFrame, setType=GOCollection())

save(H37Rv.gsc, file="H37Rv.gsc.RData")


# Write out file format for GOEAST
# Format is:
# probeID  GOIDs	annotation(optional)
# probe1	GO:0003985 // GO:0003987 // GO:0004262	annotation of probe1

#Module membership
mod_members.raw <- read.table("../../PMN_output/4.17.30_mods_members.txt",
                              head=T, sep='\t')
dim(mod_members.raw)
#[1] 3458    2

#Filter out members not in the universe:
mod_members <- mod_members.raw[mod_members.raw$gene %in% universe,]
dim(mod_members)
#[1] 2075    2, so ~1400 have no associated GO term

mod_parents <- read.table("../../PMN_output/4.17.30_mods_parsed.txt",
                          head=T, sep='\t')

#Set the universe to be all genes present in the analysis that have GO terms
my.universe <- mod_members$gene

#Only analyze modules with at least 1 meaningful probability

#Load in module statistics
mod.stats <- read.table("../../PMN_output/4.17_30mods_genes_pathsizes.txt",
                        head=T, sep='\t')
#Filter on number of probabilities > 0.4
mod.good <- mod.stats[mod.stats$thresh.0.2 > 0 & mod.stats$n_genes < 100,]$moduleID
length(mod.good)
#[1] 15 modules passing QC

mod_members.good <- mod_members[mod_members$moduleID %in% mod.good,]
dim(mod_members.good)
#[1] 389   2 561 genes in the good modules

#Create list of GO terms mapped to each gene
goeast <- file("GOEAST/GOEAST.annot.txt", "w")
write("probeID\tGOIDs", goeast, append=T)

for (gene in my.universe){
  gene.GO <- go.mtb.mappings[go.mtb.mappings$gene_id==gene,1]
  go.format <- paste(gene.GO, collapse=" // ")
  line <- paste(gene, go.format, sep="\t")
  write(line, goeast, append=T)
  
}
close(goeast)

# Write out all the modules
mods_GOEAST(mod_members.good, mod_parents)

# Regulator enrichment
#Pull out the parents. Note they are stored as a space delimited
#string (if there are 2 parents). Need to reformat into a 
#character vector of individual Rv IDs
parents.raw <- mod_parents[mod_parents$moduleID %in% mod.good,]$parents

parents <- unique(unlist(
  sapply(parents.raw, function(x) unlist(strsplit(x, " ", fixed=T)))
  ))

length(parents)
#[1] 20
write.table(parents, "GOEAST/regulators.txt", quote=F, row.names=F, col.names=F)

reg.universe <- read.table("../regulators/final.regulators")[,1]
#Create list of GO terms mapped to each gene
goeast <- file("GOEAST/regulators.annot.txt", "w")
write("probeID\tGOIDs", goeast, append=T)

for (gene in reg.universe){
  gene.GO <- go.mtb.mappings[go.mtb.mappings$gene_id==gene,1]
  go.format <- paste(gene.GO, collapse=" // ")
  line <- paste(gene, go.format, sep="\t")
  write(line, goeast, append=T)
  
}
close(goeast)



#Now do WGCNA modules
load("universe.RData")

#Load the WGCNA modules
load("../exprs/filt_pt5.net.RData")

wgcna.mod_members.raw <- data.frame(moduleID=filt_pt5.net@mergedColors,
                              gene=filt_pt5.net@peptides)
dim(wgcna.mod_members.raw)
#[1] 2158    2


#Filter out members not in the universe:
wgcna.mod_members <- wgcna.mod_members.raw[wgcna.mod_members.raw$gene %in% universe,]
dim(wgcna.mod_members)
#[1] 1336    2, so ~800 have no associated GO term

#Set the universe to be all genes present in the analysis 
#that have GO terms
wgcna.universe <- wgcna.mod_members$gene

#Create list of GO terms mapped to each gene
goeast <- file("GOEAST/WGCNA.annot.txt", "w")
write("probeID\tGOIDs", goeast, append=T)

for (gene in wgcna.universe){
  gene.GO <- go.mtb.mappings[go.mtb.mappings$gene_id==gene,1]
  go.format <- paste(gene.GO, collapse=" // ")
  line <- paste(gene, go.format, sep="\t")
  write(line, goeast, append=T)
  
}

close(goeast)


#Filter modules with the permutation test results
perm.test <- filt_pt5.net@permtest
good.mods <- perm.test[(perm.test[,5] < 0.005 & perm.test[,2] < 200),1]

wgcna.mod_members.good <- wgcna.mod_members[wgcna.mod_members$moduleID%in%good.mods,]

for (mod in unique(wgcna.mod_members.good$moduleID)){
  mod.genes <- wgcna.mod_members.good[wgcna.mod_members.good$moduleID==mod,]$gene
  fname = paste(c("WGCNA_", mod, ".txt"), collapse="")
  write.table(mod.genes, fname, col.names=F, row.names=F, quote=F)
}







