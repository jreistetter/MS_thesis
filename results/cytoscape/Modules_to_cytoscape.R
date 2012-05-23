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
  write("source\tinteraction\ttarget\tr2", out.f, append=T)
  for (i in 1:nrow(mod_adj)){
    row <- mod_adj[i,]
    r2.format <- round(as.numeric(row$r2), digits=3)
    line <- paste(c(row$gene1, row$interaction, row$gene2, r2.format), collapse="\t")
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

add_TFs <- function(adj_list, pDNA){
  adj_list$interaction <- "pp"
  tfs <- pDNA$TF
  
  for (i in 1:nrow(adj_list)){
    row <- adj_list[i,]
    
    #Check if gene1 is TF
    if (row$gene1 %in% tfs){
      tf.targets <- pDNA[pDNA$TF==row$gene1,]$target
      
      #if TF, check if the second gene is a target.
      #if so, set the appropriate interactions
      if (row$gene2 %in% tf.targets){
        adj_list[i,]$interaction <- "pd"
      }
    }
    
    #check if gene2 is TF
    if(row$gene2 %in% tfs){
      tf.targets <- pDNA[pDNA$TF==row$gene2,]$target
      
      #if TF check if second gene is target
      #if so, switch the source and target and add interaction
      if(row$gene1 %in% tf.targets){
        adj_list[i,]$gene1 <- row$gene2
        adj_list[i,]$gene2 <- row$gene1
        adj_list[i,]$interaction <- "pd"
      }
    }
  }
  return(adj_list)
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

#TF binding
tf.binding <- read.table("PMN_data/4.12.2012/H37Rv.pdna.list",
                         head=F, sep="\t")

colnames(tf.binding) <- c("TF", "target", "weight", "direction")


#Write module out to cytoscape files
mod5.adj <- mod_adj_to_list("5", wgcna.modules, filt_pt5.net, 0)

mod5.adj.tfs <- add_TFs(mod5.adj, tf.binding)
mod_adj_to_cyto(mod5.adj.tfs, "data/results/cytoscape/w.mod5.tfs.sif")


mod5.adj.annot <- annotate_adj(mod5.adj, genes.annot)


#Make node annotation file with gene name, Rv ID, PMN membership, and TF binding
w.node.annot <- merge(wgcna.modules, genes.annot,
                        by.x="gene", by.y="rvID",
                        all.x=T)

#some genes have 2 names, remove those rows
w.node.annot <- w.node.annot[-c(891, 932),]

#Add PMN module membership, not parents though
colnames(w.node.annot)[2] <- "W.module"
pmn.module.members <- pmn.modules[pmn.modules$parent==F,]

w.node.annot <- merge(w.node.annot, pmn.module.members[,c(1,3)],
                      by.x="gene", by.y="rvID",
                      all.x=T)

pmn.parents <- pmn.modules[pmn.modules$parent==T,]

w.node.annot$ParentOf <- NA
for (parent in unique(pmn.parents$rvID)){
  parent.mods.v <- pmn.parents[pmn.parents$rvID == parent,1]
  parent.mods  <- paste(parent.mods.v, collapse="|")
  if (parent %in% w.node.annot$gene){
    w.node.annot[w.node.annot$gene == parent,]$ParentOf <- parent.mods
  }
}




p.mod_membership_cyto(mod5.adj.annot, pmn.modules, 
                      "data/results/cytoscape/w.mod5.pmn.noa")

