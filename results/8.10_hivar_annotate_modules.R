#Script to map the RvIDs to common names for the module membership and parents and pathways

#Written by Joe Reistetter


library(biomaRt)
options(stringsAsFactors=F)

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
setwd("~/Dropbox/thesis_work/PMN_output/")


##############################
#
#   Annotate the modules
#
##############################
    
modules.raw <- read.table("8.10.12_30_mods_0.1_highvar_mods_members.txt",
                          head=T, sep="\t")

bac <- useMart('bacteria_mart_14', dataset='myc_30_gene')

gene.ids.all <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                  filters="tuberculist", 
                  values=modules.raw$gene, mart=bac)

nrow(modules.raw)
nrow(gene.ids.all)
# > nrow(modules.raw)
# [1] 2158
# > nrow(gene.ids.all)
# [1] 2203
# The mart query returned more rows, investigate.

sum(duplicated(gene.ids.all$tuberculist))
#[1] 48, so have dupes

# > gene.ids.all[gene.ids.all$tuberculist=="Rv3186",]
# tuberculist external_gene_id
# 74        Rv3186          Rv1370c
# 99        Rv3186           Rv2105
# 494       Rv3186          Rv3381c
# I looked these up on tbdb, all are called "transposase" but thats not their
# official symbol. Some DB mapping error, need to remove these

# Decided to use tbdb.org online annotation tool to get these:
# http://genome.tbdb.org/annotation/genome/tbdb/
#   BatchSelect.html?target=AnnotationList.html

# Write out list of Rv Ids for submission

write.table(modules.raw$gene, "../../data/exprs/high_var_gene_list.txt",
            quote=F, row.names=F, col.names=F)

gene.ids.all$tuberculist <- toupper(gene.ids.all$tuberculist)

colnames(gene.ids.all) <- c("rvID", "name")

write.table(gene.ids.all, "data/results/genes.annotated",
            col.names=T, row.names=F, sep="\t", quote=F)


good.modules <- read.table("data/results/PMN_good_modules.txt",
                           head=T, sep="\t")

colnames(good.modules) <- "moduleID"

module_members <- modules.raw[modules.raw$moduleID%in% good.modules$moduleID,]
module_members$parent <- FALSE

#Extract parents and make a new DF with parents included.
parents.raw <- read.table("PMN_output/4.17.30_mods_parsed.txt",
                          head=T, sep="\t")

parents.raw <- parents.raw[parents.raw$moduleID %in% good.modules$moduleID,]
parents <- get_parents(parents.raw)
parents$parent <- TRUE

modules <- rbind(module_members, parents)

gene.ids <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                    filters="tuberculist", 
                    values=modules$gene, mart=bac)

gene.ids$tuberculist <- toupper(gene.ids$tuberculist)
modules.annotated <- merge(modules, gene.ids, 
                           by.x="gene", 
                           by.y="tuberculist", 
                           all.x=T)[,c(2:4)]

modules.annotated <- merge(modules, gene.ids, 
                           by.x="gene", 
                           by.y="tuberculist", 
                           all.x=T)

colnames(modules.annotated) <- c("rvID", "moduleID", "parent", "name")

#Sort and change the column order
modules.annotated <- modules.annotated[with(modules.annotated, order(moduleID, -parent)),]
modules.annotated <- modules.annotated[,c(2,4,1,3)]

write.table(modules.annotated, "data/results/PMN_modules_annotated.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

##############################
#
#   Annotate the pathways
#
##############################

path.raw <- read.table("PMN_output/4.17.30_mods_pathways_for_annotation.txt",
                       head=F, sep=";")

path.raw <- path.raw[path.raw[,1] %in% good.modules$moduleID,]

path.allgenes <- unique(unlist(apply(path.raw, 1, function(x) strsplit(x[2], "\t", fixed=T))))

genes.annot <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                     filters="tuberculist", 
                     values=path.allgenes, mart=bac)

genes.annot[,1] <- toupper(genes.annot[,1])

paths.annot <- apply(path.raw, 1, function(row){
  row.genes <- unlist(strsplit(row[2], "\t", fixed=T))
  print(row.genes)
  row.annotated <- vector(mode="character", length=length(row.genes))
  for (i in c(1:length(row.genes))){
    if(nrow(genes.annot[genes.annot[,1]==row.genes[i],])==0){
      row.annotated[i] <- row.genes[i]
    }
        
    else{
    gene.name <- genes.annot[genes.annot[,1]==row.genes[i],2]
    row.annotated[i] <- gene.name
    }
  }
  return(row.annotated)
}
)


pathways <- data.frame(moduleID=path.raw[,1])
pathways$pathway <- ""

for (i in c(1:nrow(pathways))){
  pathways[i,2] <- paste(paths.annot[[i]], collapse=" -> ")
}

pathways <- pathways[with(pathways, order(moduleID)),]

write.table(pathways, "data/results/PMN_module_pathways.txt",
            row.names=F, col.names=T, sep="\t", quote=F)
