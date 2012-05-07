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
    mod.p <- hotelling.test(.~groups, data=mod.expr.groups, shrinkage=shrink, 
                            perm=perm, B=10000)
    results <- rbind(results, c(mod, mod.out$pval, mod.out$stats$statistic, 
                                mod.out$stats$nx, mod.out$stats$p))
  }
  
  return(results)
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


##ICSNP function
icsnp.p <- module.ICSNP(good.modules, modules, expr.immune, celltype)
#Throws error because of singular covariance matrix:
#Error in solve.default(S.pooled) : 
#  system is computationally singular: reciprocal condition number = 3.34324e-20

small.mods <- c("mod15", "mod16")
icsnp.p <- module.ICSNP(small.mods, modules, expr.immune, celltype)
#Works:
# mod15       mod16 
# 0.001314223 0.145787981 

#Hotelling function with S3 method
dc_mac.p <- module.Hotelling(good.modules, modules, expr.immune, celltype)
# Error in solve.default(sPooled) : 
#   system is computationally singular: reciprocal condition number = 0

dc_map.p.s3.shrink <- module.Hotelling(good.modules, modules, 
                                                   expr.immune, celltype,
                                                   shrink=T)
#Throws error:
# Error in if (denominator == 0) lambda.var = 1 else lambda.var = min(1,  : 
#   missing value where TRUE/FALSE needed

#Try using formula instead:
dc_mac.p <- module.Hotelling.formula(good.modules, modules, expr.immune, 
                                celltype, shrink=F)
# Fails correctly:
# Error in hotelling.stat(x, y, shrinkage) : 
#   The sample sizes (nx + ny) must be 1 greater than the number of columns

#Try the small modules that should work with classic Hotellings T2
dc_mac.p <- module.Hotelling.formula(small.mods, modules, expr.immune, 
                                        celltype, shrink=F)
#Matches ICSNP:
# mod15       mod16 
# 0.001314223 0.145787981 

#Run with no permutations
dc_mac.p.shrink <- module.Hotelling.formula(good.modules, modules, expr.immune, 
                                            celltype, shrink=T, perm=F)

#Gives warning because df2 is a negative number because m < p
#and df2 in Hotellings is m - p + 1
#10 warnings in all like this:
#Warning messages:
# 1: In pf(q, df1, df2, lower.tail, log.p) : NaNs produced

dc_mac.p.shrink.perm <- module.Hotelling.formula(good.modules, modules, expr.immune, 
                                celltype, shrink=T)
#Seems to work:
# mod12 mod15 mod16 mod18 mod19  mod2 mod20 mod21 mod23 mod25 mod27 mod28  mod4  mod5  mod8 
# 0.000 0.002 0.068 0.000 0.042 0.000 0.003 0.007 0.000 0.000 0.125 0.000 0.000 0.000 0.023

sum(dc_mac.p.shrink.perm==0)
#8, so not all 10 warnings produced 0 permutation p-values

dc_mac.p.shrink.adj <- p.adjust(dc_mac.p.shrink, method="BH")


#Although for example mod16 has a much lower p-value with shrinkage and permutation
#testing than the parametric test: 0.146 (parametric) vs 0.068
#mod15 on the other hand goes from 0.0013 (parametric) to 0.002


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

time.1.18.pvals <- module.Hotelling.formula(good.modules, modules, 
                                            expr.time.1.18, time.1.18$time,
                                            shrink=T)
time.1.18.pvals.adj <- p.adjust(time.1.18.pvals, method="BH")

#Try the non-shrinkage version
time.1.18.Hotelling <- module.Hotelling.formula(small.mods, modules, 
                                            expr.time.1.18, time.1.18$time,
                                            shrink=F, perm=F)
#Very small p:
# mod15        mod16 
# 3.932020e-04 5.245876e-05 

#Double check with ICSNP
time.1.18.ICSNP <- module.ICSNP(small.mods, modules, 
                                            expr.time.1.18, time.1.18$time)
#Very low pvalues:
# mod15        mod16 
# 3.932020e-04 5.245876e-05 

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

mac.aer.pvals <- module.Hotelling.formula(good.modules, modules, 
                                          expr.mac.aer, mac.aer.labels,
                                          shrink=T)
# mod12 mod15 mod16 mod18 mod19  mod2 mod20 mod21 mod23 mod25 mod27 mod28  mod4  mod5  mod8 
# 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 

mac.aer.pvals.adj <- p.adjust(mac.aer.pvals, method="BH")


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



                     
