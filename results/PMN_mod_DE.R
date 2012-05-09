# Script to perform Hotellings T2 on module members between aerobic, DC, and Mac Mtbs
# 


options(stringsAsFactors=F)

#Functions
get_module <- function(modID, modules){
  return(modules[modules$moduleID==modID,2])
}

selectGenes <- function(genes, arrays){
  genes.sel <- which(rownames(arrays) %in% genes)
  genes.expr <- t(arrays[genes.sel,])
  return(genes.expr)
}

#ICSNP package functions
calcHotellings <- function(mod.expr, grouping){
  result <- HotellingsT2(mod.expr ~ as.factor(grouping))
  return(result$p.value)
}

module.ICSNP <- function(modIDs, modules, arrays, grouping){
  results <- list()
  
  for (mod in modIDs){
    mod.genes <- get_module(mod, modules)
    mod.expr <- selectGenes(mod.genes, arrays)
    mod.p <- calcHotellings(mod.expr, grouping)
    results[[mod]] <- mod.p
  }
  
  return(unlist(results))
}


#Hotelling package functions
calc.Hotelling <- function(mod.expr, grouping){
  result <- Hott(mod.expr ~ as.factor(grouping), shrinkage=FALSE)
  return(result)
}

module.Hotelling <- function(modIDs, modules, arrays, grouping, shrink=F){
  results <- data.frame(moduleID=vector(mode="character"),
                        p=vector(mode="numeric"),
                        T_stat=vector(mode="numeric"),
                        n_samp=vector(mode="integer"),
                        n_vars=vector(mode="integer"))
  
  for (mod in modIDs){
    #Get the subset of the expression data for the module
    mod.genes <- get_module(mod, modules)
    mod.expr <- selectGenes(mod.genes, arrays)
    
    #Split the expression matrix into the two groups
    group1.id <- unique(grouping)[1]
    group2.id <- unique(grouping)[2]
    group1.rows <- which(grouping == group1.id)
    group2.rows <- which(grouping == group2.id)
    
    group1.expr <- mod.expr[group1.rows,]
    group2.expr <- mod.expr[group2.rows,]
    
    #Run the test and store the results
    mod.out <- hotelling.test(group1.expr, group2.expr, shrinkage=shrink)
    results <- rbind(results, c(mod, mod.out$pval, mod.out$stats$statistic, 
                                mod.out$stats$nx, mod.out$stats$p))
  }
  
  return(results)
}

module.Hotelling.formula <- function(modIDs, modules, arrays, grouping, 
                                     shrink=F, perm=T){
  results <- data.frame(moduleID=vector(mode="character"),
                        p=vector(mode="numeric"),
                        T_stat=vector(mode="numeric"),
                        n_samp=vector(mode="integer"),
                        n_vars=vector(mode="integer"))
  
  for (mod in modIDs){
    #Get the subset of the expression data for the module
    mod.genes <- get_module(mod, modules)
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
immune.time.DE <- function(samples, arrays, good.modules, modules, time){
  #Macs vs DCs at 18h
  samp.sel <- samples[(samples$celltype != "Aerobic" & samples$time==time),]
  expr <- arrays[,colnames(arrays)%in%samp.sel$filename]

  pvals <- module.Hotelling.formula(good.modules, modules, 
                                               expr, samp.sel$celltype,
                                               shrink=T)
  
  pvals$p.adj <- p.adjust(pvals$p, method="BH")
  return(pvals)
}

library(ICSNP)
library(Hotelling)

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



dc_mac.p.shrink.perm <- module.Hotelling.formula(good.modules, modules, expr.immune, 
                                celltype, shrink=T)
dc_mac.p.shrink.perm$p.adj <- p.adjust(dc_mac.p.shrink.perm$p, method="BH")
write.table(dc_mac.p.shrink.perm, "data/results/DC_vs_Mac_DE.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

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

time.1.18.pvals <- module.Hotelling.formula(good.modules, modules, 
                                            expr.time.1.18, time.1.18$time,
                                            shrink=T)
time.1.18.pvals$p.adj <- p.adjust(time.1.18.pvals$p, method="BH")
write.table(time.1.18.pvals, "data/results/Immunes_1h_vs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs vs DCs at 18h
dc_mac_18h <- BUGS58.samples[(BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time =="18h"),]
expr.dc_mac_18h <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%dc_mac_18h$filename]
dim(expr.dc_mac_18h)
#[1] 3765   12

dc_mac_18h.pvals <- module.Hotelling.formula(good.modules, modules, 
                                            expr.dc_mac_18h, dc_mac_18h$celltype,
                                            shrink=T)
dc_mac_18h.pvals$p.adj <- p.adjust(dc_mac_18h.pvals$p, method="BH")

write.table(dc_mac_18h.pvals, "data/results/DC_vs_macs_18h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

#Macs vs DCs at 4h
dc_mac_4h <- BUGS58.samples[(BUGS58.samples$celltype != "Aerobic" & BUGS58.samples$time =="4h"),]
expr.dc_mac_4h <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%dc_mac_4h$filename]
dim(expr.dc_mac_4h)
#[1] 3765   12

dc_mac_4h.pvals <- module.Hotelling.formula(good.modules, modules, 
                                            expr.dc_mac_4h, dc_mac_4h$celltype,
                                            shrink=T)
dc_mac_4h.pvals$p.adj <- p.adjust(dc_mac_4h.pvals$p, method="BH")

write.table(dc_mac_4h.pvals, "data/results/DC_vs_macs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)



dc_mac_4h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, good.modules,
                                  modules, "4h")

write.table(dc_mac_4h.pvals, "data/results/DC_vs_macs_4h.txt",
            col.names=T, sep="\t", quote=F, row.names=F)

dc_mac_1h.pvals <- immune.time.DE(BUGS58.samples, BUGS58.arrays, good.modules,
                                  modules, "1h")

#### 
#
#    MACS
#
####

snp.alleles <- data.frame(snpID=allSNPIDs)

#Macs vs Aerobic
mac.aer.arr <- BUGS58.samples[BUGS58.samples$celltype != "DC",1]
mac.aer.labels <- BUGS58.samples[BUGS58.samples$celltype != "DC", 2]
expr.mac.aer <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.aer.arr]
dim(expr.mac.aer)
#[1] 3765   26, 26 total arrays

mac.aer.pvals <- module.Hotelling.formula(good.modules, modules, 
                                          expr.mac.aer, mac.aer.labels,
                                          shrink=T)
mac.aer.pvals$p.adj <- p.adjust(mac.aer.pvals, method="BH")


#Macs 0h vs 4h
mac.1.4 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "4h"),]
expr.mac.1.4 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.4$filename]
dim(expr.mac.1.4)
#[1] 3765   12, 12 arrays

mac.1.4.pvals <- module.Hotelling.formula(good.modules, modules, 
                                          expr.mac.1.4, mac.1.4$time,
                                          shrink=T)


#Macs 0h vs 18 h
mac.1.18 <- BUGS58.samples[BUGS58.samples$celltype == "MDM" & BUGS58.samples$time %in% c("1h", "18h"),]
expr.mac.1.18 <- BUGS58.arrays[,colnames(BUGS58.arrays)%in%mac.1.18$filename]
dim(expr.mac.1.18)
#[1] 3765   12, 12 arrays

mac.1.18.pvals <- module.Hotelling.formula(good.modules, modules, 
                                           expr.mac.1.18, mac.1.18$time,
                                           shrink=T)


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

DC.aer.pvals <- moduleHotellings(good.modules, modules, 
                                 expr.DC.aer, DC.aer.labels,
                                 shrink=T)

DC.aer.pvals.adj <- p.adjust(DC.aer.pvals, method="BH")



                     
