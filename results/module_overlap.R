# Script to calculate how many members overlap between
# PMN and WGCNA modules
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

#Functions

make_c_table <- function(pmn.id, wgcna.id, pmn.modules, wgcna.modules, n){
  # n = total number of genes
  # m = number in common between n1 and n2
  # n1 = pmn.genes
  # n2 = wgcna.genes
  # 
  #    Table:
  #
  #     m         n1 - m
  #    n2 - m    n-n1-n2+m
  #
  pmn.mod <- get_module(pmn.id, pmn.modules)
  wgcna.mod <- get_module(wgcna.id, wgcna.modules)
  m <- length(intersect(pmn.mod, wgcna.mod))
  n1 <- sum(pmn.mod %in% wgcna.modules$gene)
  n2 <- length(wgcna.mod)
  
  c_table <- matrix(c(m, n1 - m, n2 - m, n - n1 - n2 + m),
                    nrow=2, byrow=T)
  
  return(c_table)
}

calc_overlap <- function(pmn.id, wgcna.id, pmn.modules, wgcna.modules, n){
  c_table <- make_c_table(pmn.id, wgcna.id, pmn.modules, wgcna.modules, n)
  p.val <- fisher.test(c_table)$p.value
  
  return(p.val)
  
}

overlap_table <- function(pmn.mod.ids, wgcna.mod.ids, pmn.modules, wgcna.modules, n){
  pTable <- matrix(0, nrow=length(pmn.mod.ids), ncol=length(wgcna.mod.ids))
  
  for (i in 1:length(pmn.mod.ids)){
    pmn.id <- pmn.mod.ids[i]
    for (j in 1:length(wgcna.mod.ids)){
      wgcna.id <- wgcna.mod.ids[j]
      overlap.p <- calc_overlap(pmn.id, wgcna.id, pmn.modules, wgcna.modules, n)
      pTable[i, j] <- overlap.p
    }
  }
  
  p.adj <- p.adjust(c(pTable), method="BH")
  pTable.adj <- matrix(p.adj, ncol=ncol(pTable))
  
  rownames(pTable.adj) <- pmn.mod.ids
  colnames(pTable.adj) <- wgcna.mod.ids
  
  return(-log10(pTable.adj))
}



get_module <- function(modID, modules){
  return(modules[modules$moduleID==modID,2])
}


#########################
#
#  Main
#
#########################

library(gplots)

#Load data

setwd("~/Dropbox/thesis_work/data/")

#Load PMN modules
pmn.modules <- read.table("../PMN_output/4.17.30_mods_members.txt",
                  head=T, sep='\t')

pmn.mod.stats <- read.table("../PMN_output/4.17_30mods_genes_pathsizes.txt",
                            head=T, sep='\t')

pmn.modIDs <- pmn.mod.stats[pmn.mod.stats$n_genes < 100 & pmn.mod.stats$thresh.0.2 > 0,]$moduleID

write.table(pmn.modIDs, 
            "results/PMN_good_modules.txt", sep="\t", 
            row.names=F, col.names=T, quote=F)


#Load WGCNA modules
load("exprs/filt_pt5.net.RData")
wgcna.stats <- as.data.frame(filt_pt5.net@permtest)
wgcna.stats$p.adj <- p.adjust(wgcna.stats[,5], method="BH")
write.table(wgcna.stats, "results/WGCNA_module_stats.txt",
            row.names=F, col.names=T, sep="\t", quote=F)
wgcna.modules <- data.frame(moduleID=filt_pt5.net@mergedColors, gene=filt_pt5.net@peptides)
wgcna.modIDs <- wgcna.stats[wgcna.stats[,2] < 200 & wgcna.stats[,6] < 0.05,1]
write.table(wgcna.modIDs, "results/WGCNA_good_modules.txt",
            row.names=F, col.names=T, quote=F, sep="\t")




#Calculate overlap
n <- length(unique(c(pmn.modules$gene, wgcna.modules$gene)))

overlap.table <- overlap_table(pmn.modIDs, wgcna.modIDs, pmn.modules, wgcna.modules, n)

write.table(overlap.table, "results/Module_overlap.txt",
            row.names=T,
            col.names=T,
            quote=F,
            sep="\t")

overlap.heat <- apply(overlap.table, c(1,2), function(x){
  if (x < 2){
    return(0)
  }
  return(x)
}
                      )
                      

heatmap.2(overlap.heat,
          col=colorRampPalette(c("white", "darkblue"))(100),
          trace="none",
          Rowv=F,
          Colv=F,
          scale="none",
          xlab="WGCNA modules",
          ylab="PMN modules",
          density.info="none")

