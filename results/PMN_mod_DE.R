# Script to perform Hotellings T2 on module members between aerobic, DC, and Mac Mtbs
# and make heatmaps of module expression


options(stringsAsFactors=F)

#Functions
get_module <- function(modID, modules){
  return(modules[modules$moduleID==modID,2])
}

get_parents <- function(modID, parents){
  raw <- parents[parents$moduleID==modID,2]
  parent.genes <- unique(unlist(strsplit(raw, " ", fixed=T)))
  return(parent.genes)
}

selectGenes <- function(genes, arrays){
  genes.sel <- which(rownames(arrays) %in% genes)
  genes.expr <- t(arrays[genes.sel,])
  return(genes.expr)
}

#Hotelling package functions
calc.Hotelling <- function(mod.expr, grouping){
  result <- Hott(mod.expr ~ as.factor(grouping), shrinkage=FALSE)
  return(result)
}


module.Hotelling <- function(modIDs, modules, parents, arrays, grouping, 
                                     shrink=F, perm=T){
  results <- data.frame(moduleID=vector(mode="character"),
                        p=vector(mode="numeric"),
                        T_stat=vector(mode="numeric"),
                        n_samp=vector(mode="integer"),
                        n_vars=vector(mode="integer"))
  
  for (mod in modIDs){
    #Get the subset of the expression data for the module
    mod.genes <- get_module(mod, modules)
    mod.genes <- c(mod.genes, get_parents(mod, parents))
    mod.expr <- as.data.frame(selectGenes(mod.genes, arrays))
    
    #Create group vector of {1,2}
    groups <- rep(1, length(grouping))
    groups[grouping==unique(grouping)[2]] <- 2
    
    #Add to matrix
    mod.expr.groups <- cbind(groups, mod.expr)
    
    #Run the test and store the results
    mod.out <- hotelling.test(.~groups, data=mod.expr.groups, shrinkage=shrink, 
                            perm=perm, B=10000)
    results <- rbind(results, c(mod, mod.out$pval, mod.out$stats$statistic, 
                                mod.out$stats$nx, mod.out$stats$p))
  }
  colnames(results) <- c("moduleID", "p", "T_stat", "n_samp", "n_vars")
  return(results)
}

#Function to do the immune cells against each other at each timepoint
immune.time.DE <- function(samples, arrays, good.modules, modules, parents, time){
  #Macs vs DCs at 18h
  samp.sel <- samples[(samples$celltype != "Aerobic" & samples$time==time),]
  expr <- arrays[,colnames(arrays)%in%samp.sel$filename]

  pvals <- module.Hotelling(good.modules, modules, parents,
                                               expr, samp.sel$celltype,
                                               shrink=T)
  
  pvals$p.adj <- p.adjust(pvals$p, method="BH")
  return(pvals)
}

heat_labels <- function(arrayIDs, samples){
  # arrayIDs - char vector of filenames of arrays
  # samples - dataframe of the sample metadata
  
  array.info <- samples[samples$filename %in% arrayIDs,]
  array.info$labels <- unlist(apply(array.info, 1, function(array){
    label <- paste(c(array[2],
                     array[3],
                     array[4]),
                   collapse=" - "
                   )
    return(label)
  }
                                    ))
  return(array.info$labels)
}

mod.heat <- function(mod.genes, expr, samples, title){
  # mod.genes - char vector of genes for heatmap
  # expr - all expression data
  mod.expr <- selectGenes(mod.genes, expr)
  #mod.expr <- mod.expr[complete.cases(mod.expr),]
  mod.expr.t <- t(mod.expr)
  colnames(mod.expr.t) <- heat_labels(rownames(mod.expr), samples)
  #Set colors for DCs and pass to function.
  cell.type <- grepl("DC", colnames(mod.expr.t), fixed=T)
  col.side <- rep("seagreen", length(cell.type))
  col.side[cell.type] <- "skyblue"
  print(col.side)
  heatmap.2(mod.expr.t,
            cexRow = 0.5,
            cexCol = 0.5,
            na.rm=T, 
            trace="none",
            symkey=T,
            col=redgreen,
            key=TRUE,
            ColSideColors=col.side)
}

mod.heat.time <- function(mod.genes, expr, samples, title){
  # mod.genes - char vector of genes for heatmap
  # expr - all expression data
  mod.expr <- selectGenes(mod.genes, expr)
  #mod.expr <- mod.expr[complete.cases(mod.expr),]
  mod.expr.t <- t(mod.expr)
  colnames(mod.expr.t) <- heat_labels(rownames(mod.expr), samples)
  print(colnames(mod.expr.t))
  #Set colors for DCs and pass to function.
  time.4h <- grepl("4h", colnames(mod.expr.t), fixed=T)
  time.18h <- grepl("18h", colnames(mod.expr.t), fixed=T)
  col.side <- rep("orange", length(time.4h))
  col.side[time.4h] <- "purple"
  col.side[time.18h] <- "red"
  print(col.side)
  heatmap.2(mod.expr.t,
            cexRow = 0.5,
            cexCol = 0.5,
            na.rm=T, 
            trace="none",
            symkey=T,
            col=redgreen,
            key=TRUE,
            ColSideColors=col.side)
}

library(gplots)
library(Hotelling)

setwd("~/Dropbox/thesis_work/")

#Load data

#Module data

#Filter out the modules with no probabilities
modules <- read.table("PMN_output/4.17.30_mods_members.txt",
                          head=T, sep="\t")

parents <- read.table("PMN_output/4.17.30_mods_parsed.txt",
                      head=T, sep="\t")

modules.stats <- read.table("PMN_output/4.17_30mods_genes_pathsizes.txt",
                            head=T, sep="\t")

good.modules <- modules.stats[modules.stats$thresh.0.2 > 0 & modules.stats$n_genes < 100,]$moduleID

#Arrays
load("data/exprs/EBUGS58/EBUGS58.arrays.RData")
#Samples
load("data/exprs/EBUGS58/BUGS58.samples.RData")



##############################
#
#     Hotellings T2 Analysis
#
##############################

## DC vs Macs
immune.arrays <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic",1]
celltype <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic", 2]
expr.immune <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%immune.arrays]
dim(expr.immune)
#[1] 3765   36, 36 total arrays

dc_mac.p.shrink.perm <- module.Hotelling(good.modules, modules, parents, expr.immune, 
                                celltype, shrink=T)
dc_mac.p.shrink.perm$p.adj <- p.adjust(dc_mac.p.shrink.perm$p, method="BH")
write.table(dc_mac.p.shrink.perm, "data/results/PMN_DC_vs_Mac_DE.txt",
            col.names=T, sep="\t", quote=F, row.names=F)





mod2 <- get_module("mod2", modules)
mod.heat(mod2, expr.immune, BUGS58.samples)
## DCs = blue, Macs = green
mod.heat.time(mod2, expr.immune, BUGS58.samples)
## 1h = orange, 4h = purple, 18h = red




colnames(e.t) <- heat_labels(colnames(expr.immune), BUGS58.samples)
heatmap.2(e.t)



#### 
#
#   Time
#
####

# 1h vs 18h, Immune cells
time.1.18 <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.time.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%time.1.18$filename]
dim(expr.time.1.18)
#[1] 3765   24, 24 arrays

time.1.18.pvals <- module.Hotelling(good.modules, modules, parents,
                                            expr.time.1.18, time.1.18$time,
                                            shrink=T)
time.1.18.pvals$p.adj <- p.adjust(time.1.18.pvals$p, method="BH")
write.table(time.1.18.pvals, "data/results/PMN_immunes_1h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs vs DCs at 18h
dc_mac_18h <- BUGS58.samples[(BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time =="18h"),]
expr.dc_mac_18h <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%dc_mac_18h$filename]
dim(expr.dc_mac_18h)
#[1] 3765   12

dc_mac_18h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, good.modules, 
                                   modules, parents, "18h")

write.table(dc_mac_18h.pvals, "data/results/PMN_DC_vs_macs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

dc_mac_4h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, good.modules,
                                  modules, parents, "4h")

write.table(dc_mac_4h.pvals, "data/results/PMN_DC_vs_macs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

dc_mac_1h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, good.modules,
                                  modules, parents, "1h")

mods.1h <- good.modules[-which(good.modules %in% c("mod18"))]

dc_mac_1h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, mods.1h,
                                  modules, parents, "1h")

write.table(dc_mac_1h.pvals, "data/results/PMN_DC_vs_macs_1h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)


#Merge into one big DF
colnames(dc_mac_1h.pvals)[6] <- "p.1h"
colnames(dc_mac_4h.pvals)[6] <- "p.4h"
colnames(dc_mac_18h.pvals)[6] <- "p.18h"

dc_mac_time_p <- merge(dc_mac_1h.pvals, dc_mac_4h.pvals)

#### 
#
#    MACS
#
####

#Macs vs Aerobic
mac.aer.arr <- BUGS58.samples[BUGS58.samples$celltype != "DC",1]
mac.aer.labels <- BUGS58.samples[BUGS58.samples$celltype != "DC", 2]
expr.mac.aer <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.aer.arr]
dim(expr.mac.aer)
#[1] 3765   26, 26 total arrays

mac.aer.pvals <- module.Hotelling(good.modules, modules, parents,
                                          expr.mac.aer, mac.aer.labels,
                                          shrink=T)
mac.aer.pvals$p.adj <- p.adjust(mac.aer.pvals$p, method="BH")

write.table(mac.aer.pvals, "data/results/PMN_Macs_vs_Aerobic.txt",
            col.names=T, sep="\t", quote=F, row.names=F)


#Macs 1h vs 4h
mods.1.4 <- good.modules[-which(good.modules %in% c("mod18"))]


mac.1.4 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "4h"),]
expr.mac.1.4 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.4$filename]
dim(expr.mac.1.4)
#[1] 3765   12, 12 arrays

mac.1.4.pvals <- module.Hotelling(mods.1.4, modules, parents,
                                          expr.mac.1.4, mac.1.4$time,
                                          shrink=T)

mac.1.4.pvals$p.adj <- p.adjust(mac.1.4.pvals$p, method="BH")

write.table(mac.1.4.pvals,"data/results/PMN_Macs_1h_vs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs 4h vs 18 h
mac.4.18 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("4h", "18h"),]
expr.mac.4.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.4.18$filename]
dim(expr.mac.4.18)
#[1] 3765   12, 12 arrays

mac.4.18.pvals <- module.Hotelling(mods.1.4, modules, parents,
                                   expr.mac.4.18, mac.4.18$time,
                                   shrink=T)

mac.4.18.pvals$p.adj <- p.adjust(mac.4.18.pvals$p, method="BH")

write.table(mac.4.18.pvals,"data/results/PMN_Macs_4h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs 1h vs 18 h
mac.1.18 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.mac.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.18$filename]
dim(expr.mac.1.18)
#[1] 3765   12, 12 arrays

mac.1.18.pvals <- module.Hotelling(mods.1.4, modules, parents,
                                           expr.mac.1.18, mac.1.18$time,
                                           shrink=T)

mac.1.18.pvals$p.adj <- p.adjust(mac.1.18.pvals$p, method="BH")

write.table(mac.1.18.pvals,"data/results/PMN_Macs_1h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#### 
#
#    DCs
#
####
#DC vs aerobic
DC.aer.arr <- BUGS58.samples[BUGS58.samples$celltype != "MDM",1]
DC.aer.labels <- BUGS58.samples[BUGS58.samples$celltype != "MDM", 2]
expr.DC.aer <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%DC.aer.arr]
dim(expr.DC.aer)
#[1] 3765   26, 26 total arrays

DC.aer.pvals <- module.Hotelling(good.modules, modules, parents,
                                 expr.DC.aer, DC.aer.labels,
                                 shrink=T)

DC.aer.pvals$adj <- p.adjust(DC.aer.pvals$p, method="BH")

write.table(DC.aer.pvals,"data/results/PMN_DCs_vs_Aerobic.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#DCs 1h vs 4h
mods.1.4 <- good.modules[-which(good.modules %in% c("mod18"))]


DC.1.4 <- BUGS58.samples[BUGS58.samples$celltype == "DC" & BUGS58.samples$time %in% c("1h", "4h"),]
expr.DC.1.4 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%DC.1.4$filename]
dim(expr.DC.1.4)
#[1] 3765   12, 12 arrays

DC.1.4.pvals <- module.Hotelling(good.modules, modules, parents,
                                  expr.DC.1.4, DC.1.4$time,
                                  shrink=T)

DC.1.4.pvals$p.adj <- p.adjust(DC.1.4.pvals$p, method="BH")

write.table(DC.1.4.pvals,"data/results/PMN_DCs_1h_vs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)


#DC 4h vs 18h
DC.4.18 <- BUGS58.samples[BUGS58.samples$celltype == "DC" & BUGS58.samples$time %in% c("4h", "18h"),]
expr.DC.4.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%DC.4.18$filename]
dim(expr.DC.4.18)
#[1] 3765   12, 12 arrays

DC.4.18.pvals <- module.Hotelling(good.modules, modules, parents,
                                 expr.DC.4.18, DC.4.18$time,
                                 shrink=T)

DC.4.18.pvals$p.adj <- p.adjust(DC.4.18.pvals$p, method="BH")

write.table(DC.4.18.pvals,"data/results/PMN_DCs_4h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#DC 1h vs 18h
DC.1.18 <- BUGS58.samples[BUGS58.samples$celltype == "DC" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.DC.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%DC.1.18$filename]
dim(expr.DC.1.18)
#[1] 3765   12, 12 arrays

DC.1.18.pvals <- module.Hotelling(good.modules, modules, parents,
                                  expr.DC.1.18, DC.1.18$time,
                                  shrink=T)

DC.1.18.pvals$p.adj <- p.adjust(DC.1.18.pvals$p, method="BH")

write.table(DC.1.18.pvals,"data/results/PMN_DCs_1h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)





#Old code
dc_mac_18h.pvals$p.adj <- p.adjust(dc_mac_18h.pvals$p, method="BH")

write.table(dc_mac_18h.pvals, "data/results/DC_vs_macs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs vs DCs at 4h
dc_mac_4h <- BUGS58.samples[(BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time =="4h"),]
expr.dc_mac_4h <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%dc_mac_4h$filename]
dim(expr.dc_mac_4h)
#[1] 3765   12

dc_mac_4h.pvals <- module.Hotelling(good.modules, modules, parents,
                                    expr.dc_mac_4h, dc_mac_4h$celltype,
                                    shrink=T)
dc_mac_4h.pvals$p.adj <- p.adjust(dc_mac_4h.pvals$p, method="BH")

write.table(dc_mac_4h.pvals, "data/results/DC_vs_macs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)


                     
