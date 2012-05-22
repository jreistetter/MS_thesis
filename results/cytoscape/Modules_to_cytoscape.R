# Script to export the WGCNA modules to cytoscape format.
# For each node, will include:
#       -PMN module membership
#       -TFs that regulate the node
#       -If the node is a parent of a module

#Functions
setwd("~/Dropbox/thesis_work/")
source("./code/thesis_funcs.R")

get_gene_adj <- function(gene1, gene2, adj){
  # Retrieves r^2 from WGCNA adjacency matrix
  gene1.idx <- which(rownames(adj)==gene1)
  gene2.idx <- which(colnames(adj)==gene2)
  adj.genes <- adj[gene1.idx, gene2.idx]
  return(adj.genes)
}


mod_adj_to_list <- function(modID, modules, net_obj, threshold){
  adj <- net_obj@adjacency
  mod <- get_module(modID, modules)
  out.df <- data.frame(gene1=c(), gene2=c(), r2=c())
  
  while (length(mod) > 1){
    gene1 <- mod[1]
    for (gene in mod[2:length(mod)]){
      gene.adj <- get_gene_adj(gene1, gene, adj)
      if (gene.adj > threshold){
        out.df <- rbind(out.df, c(gene1, gene, gene.adj))
      }
    }
    
    mod <- mod[2:length(mod)]
  }
  
  colnames(out.df) <- c("gene1", "gene2", "r2")
  return(out.df)
  
}

mod_adj_to_cyto <- function(mod_adj, out_path){
  out.f <- file(out_path, "w")
  write("source\ttarget\tr2\trv1\trv2", out.f, append=T)
  for (i in 1:nrow(mod_adj)){
    row <- mod_adj[i,]
    r2.format <- round(as.numeric(row[3]), digits=3)
    line <- paste(c(row[4], row[5], r2.format, row[1], row[2]), collapse="\t")
    write(line, out.f, append=T)
  }
  close(out.f)
}

p.mod_membership_cyto <- function(w.mod_adj, p.modules, out_path){
  w.genes <- unique(c(w.mod_adj$name1, w.mod_adj$name2))
  
  out.f <- file(out_path, "w")
  write("PMN_module", out.f, append=T)
  
  for (gene in w.genes){
    if (!(gene %in% p.modules$name)){
      line <- paste(c(gene, "NA"), collapse = " = ")
      write(line, out.f, append=T)
    }
    else{
      gene.mod <- p.modules[p.modules$name==gene & p.modules$parent==F,]$moduleID
      line <- paste(c(gene, gene.mod), collapse = " = ")
      write(line, out.f, append=T)
    }
  }
  
  close(out.f) 
}

annotate_adj <- function(adj.df, annot.df){
  adj.df$name1 <- ""
  adj.df$name2 <- ""
  for (i in 1:nrow(adj.df)){
    gene1.annot <- annot.df[annot.df$rvID == adj.df[i,]$gene1,]$name
    gene2.annot <- annot.df[annot.df$rvID == adj.df[i,]$gene2,]$name
    adj.df[i,]$name1 <- gene1.annot
    adj.df[i,]$name2 <- gene2.annot
  }
  return(adj.df)
}


# main
options(stringsAsFactors=F)

#Load data

#gene annotations
genes.annot <- read.table("data/results/genes.annotated",
                          head=T, sep="\t")
# PMN module membership
pmn.modules <- read.table("data/results/PMN_modules_annotated.txt",
                          head=T, sep='\t')

#WGCNA modules
load("data/exprs/filt_pt5.net.RData")
wgcna.modules.raw <- data.frame(moduleID=filt_pt5.net@mergedColors, 
                            gene=filt_pt5.net@peptides)

wgcna.good <- read.table("data/results/WGCNA_good_modules.txt",
                         head=T)[,1]

wgcna.modules <- wgcna.modules.raw[wgcna.modules.raw$moduleID%in%wgcna.good,]
w.adj <- filt_pt5.net@adjacency



mod5.adj <- mod_adj_to_list("5", wgcna.modules, filt_pt5.net, 0)
mod5.adj.annot <- annotate_adj(mod5.adj, genes.annot)
mod_adj_to_cyto(mod5.adj.annot, "data/results/cytoscape/w.mod5.sif")
p.mod_membership_cyto(mod5.adj.annot, pmn.modules, 
                      "data/results/cytoscape/w.mod5.pmn.noa")

