# Script to perform Hotellings T2 on module members between aerobic, DC, and Mac Mtbs
# 
# Written by Joe Reistetter

options(stringsAsFactors=F)

selectGenes <- function(genes, arrays){
  genes.sel <- which(rownames(arrays) %in% genes)
  genes.expr <- t(arrays[genes.sel,])
  return(genes.expr)
}

calcHotellings <- function(mod.expr, grouping){
  result <- HotellingsT2(mod.expr ~ as.factor(grouping))
  return(result$p.value)
}

moduleHotellings <- function(modIDs, modules, arrays, grouping){
  results <- list()
  
  for (mod in modIDs){
    mod.genes <- modules[modules$moduleID == mod, 2]
    mod.expr <- selectGenes(mod.genes, arrays)
    mod.p <- calcHotellings(mod.expr, grouping)
    results[[mod]] <- mod.p
  }
  
  return(unlist(results))
}

#Testing Hotelling package
calc2 <- function(mod.expr, grouping){
  result <- hotelling.test(mod.expr ~ as.factor(grouping))
  return(result$pval)
}

module2 <- function(modIDs, modules, arrays, grouping){
  results <- list()
  
  for (mod in modIDs){
    mod.genes <- modules[modules$moduleID == mod, 2]
    mod.expr <- selectGenes(mod.genes, arrays)
    mod.p <- calc2(mod.expr, grouping)
    results[[mod]] <- mod.p
  }
  
  return(unlist(results))
}

library(ICSNP)

setwd("~/Dropbox/thesis_work/")

#Load data

#Module data

#Filter out the modules with no probabilities
modules <- read.table("PMN_output/4.17.30_mods_members.txt",
                          head=T, sep="\t")

modules.stats <- read.table("PMN_output/4.17_30mods_genes_pathsizes.txt",
                            head=T, sep="\t")


#Arrays
load("data/exprs/EBUGS58/EBUGS58.arrays.RData")
#Samples
load("data/exprs/EBUGS58/BUGS58.samples.RData")



##############################
#
#     Hotellings T2 Analysis
#
##############################

good.modules <- modules.stats[modules.stats$thresh.0.2 > 0 & modules.stats$n_genes < 100,]$moduleID

## DC vs Macs
immune.arrays <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic",1]
celltype <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic", 2]
expr.immune <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%immune.arrays]
dim(expr.immune)
#[1] 3765   36, 36 total arrays

#Some modules throw error when test is done
immune.modules <- good.modules[-which(good.modules%in%c("mod18", "mod12", "mod28"))]

immune.pvals <- moduleHotellings(immune.modules, modules, expr.immune, celltype)
immune.pvals.adj <- p.adjust(immune.pvals, method="BH")

p <- module2(good.modules, modules, expr.immune, celltype)


#### 
#
#   Time
#
####

# 1h vs 18h
time.1.18 <- BUGS58.samples[BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.time.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%time.1.18$filename]
dim(expr.time.1.18)
#[1] 3765   24, 24 arrays

time.1.18.modules <- good.modules[-which(good.modules%in%c("mod12", "mod18", "mod20",
                                                           "mod23", "mod25", "mod28",
                                                           "mod4", "mod5"))]

time.1.18.pvals <- moduleHotellings(time.1.18.modules, modules, expr.time.1.18, time.1.18$time)
time.1.18.pvals.adj <- p.adjust(time.1.18.pvals, method="BH")

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

mac.aer.modules <- good.modules[-which(good.modules%in%c("mod12", "mod18", "mod20", 
                                                         "mod25", "mod28", "mod4"))]

mac.aer.pvals <- moduleHotellings(mac.aer.modules, modules, expr.mac.aer, mac.aer.labels)
mac.aer.pvals.adj <- p.adjust(mac.aer.pvals, method="BH")

#Macs 0h vs 4h
mac.1.4 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "4h"),]
expr.mac.1.4 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.4$filename]
dim(expr.mac.1.4)
#[1] 3765   12, 12 arrays

mac.1.4.modules <- good.modules[-which(good.modules%in%c("mod12", "mod18", "mod2",
                                                         "mod20", "mod21", "mod23",
                                                         "mod25", "mod27", "mod28",
                                                         "mod4", "mod5"))]

mac.1.4.pvals <- moduleHotellings(mac.1.4.modules, modules, expr.mac.1.4, mac.1.4$time)
t1 <- module2(good.modules, modules, expr.mac.1.4, mac.1.4$time)

#Macs 0h vs 18 h
mac.1.18 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.mac.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.18$filename]
dim(expr.mac.1.18)
#[1] 3765   12, 12 arrays

#No mods work
mac.1.18.pvals <- moduleHotellings(mac.1.18.modules, modules, expr.mac.1.18, mac.1.18$time)

#DC vs aerobic
DC.aer.arr <- BUGS58.samples[BUGS58.samples$celltype != "MDM",1]
DC.aer.labels <- BUGS58.samples[BUGS58.samples$celltype != "MDM", 2]
expr.DC.aer <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%DC.aer.arr]
dim(expr.DC.aer)
#[1] 3765   26, 26 total arrays

DC.aer.modules <- good.modules[-which(good.modules%in%c("mod12", "mod18", "mod20",
                                                        "mod25", "mod28", "mod4"))]

DC.aer.pvals <- moduleHotellings(DC.aer.modules, modules, expr.DC.aer, DC.aer.labels)
DC.aer.pvals.adj <- p.adjust(DC.aer.pvals, method="BH")



                     
