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

# main
options(stringsAsFactors=F)

#Load data

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
  write("source\ttarget\tr2", out.f, append=T)
  for (i in 1:nrow(mod_adj)){
    row <- mod_adj[i,]
    r2.format <- round(as.numeric(row[3]), digits=3)
    line <- paste(c(row[1], row[2], r2.format), collapse="\t")
    write(line, out.f, append=T)
  }
  close(out.f)
}

p.mod_membership_cyto <- function(w.mod_adj, p.modules, out_path){
 w.genes <- unique(c(w.mod_adj$gene1, w.mod_adj$gene2))
 
 out.f <- file(out_path, "w")
 write("PMN_module", out.f, append=T)
 
 for (gene in w.genes){
   if (!(gene %in% p.modules$rvID)){
     line <- paste(c(gene, "NA"), collapse = " = ")
     write(line, out.f, append=T)
   }
   else{
     gene.mod <- p.modules[p.modules$rvID==gene & p.modules$parent==F,]$moduleID
     line <- paste(c(gene, gene.mod), collapse = " = ")
     write(line, out.f, append=T)
   }
 }
 
 close(out.f) 
}

mod5.adj <- mod_adj_to_list("5", wgcna.modules, filt_pt5.net, 0)
mod_adj_to_cyto(mod5.adj, "data/results/cytoscape/w.mod5.sif")
p.mod_membership_cyto(mod5.adj, pmn.modules, 
                      "data/results/cytoscape/w.mod5.pmn.noa")


mod2 <- pmn.modules[pmn.modules$moduleID=="mod2",]$gene
w5 <- wgcna.modules[wgcna.modules$moduleID=="5",]$gene

w5.net <- data.frame(gene1=c(w5[1]), gene2=c(w5[2]))

w5.net$adj <- get_gene_adj("RV0001", "RV0002", filt_pt5.net)

for (i in 2:(length(w5)-1)){
  gene1 <- w5[i]
  for (j in i+1:length(w5)){
    gene2 <- w5[j]
    gene.adj <- get_gene_adj(gene1, gene2, filt_pt5.net)
    w5.net <- rbind(w5.net, c(gene1, gene2, gene.adj))
  }
}

w5.filt <- w5.net[w5.net$adj > ]