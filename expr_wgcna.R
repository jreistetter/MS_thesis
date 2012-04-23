# Script to build a WGCNA in expression data to see if
# expression data might be too noisy or bad for PMN
# 
# Written by Joe Reistetter

#Load data
#Laptop
setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/")

load("expr.mean.RData")

#Functions want samples as rows, genes as columns so need to transpose
expr.t <- t(expr.mean)

library(WGCNA)
pickSoftThreshold(data=expr.mean, networkType="signed")

library(procona)
net <- newProconaObj(networkName="Mtb", 
                     pepdat=expr.t, scaleFreeThreshold=0.7, signed="signed")

save(net, file="WGCNA_net1.RData")

That produces a network object with a bunch of different bits…
in the terminal type "net@"  then hit the tab key a couple of times to 
see them all.

To get the number of modules formed:
  
  table(net@mergedColors)
