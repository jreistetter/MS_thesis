# Script to build a WGCNA in expression data to see if
# expression data might be too noisy or bad for PMN
# 
# Written by Joe Reistetter

#Load data
#Laptop
#setwd("~/schoolDB/Dropbox/thesis_work/data/exprs/")
#OHSU
setwd("~/Dropbox/thesis_work/data/exprs/")


load("expr.mean.RData")

#Functions want samples as rows, genes as columns so need to transpose
expr.t <- t(expr.mean)

library(WGCNA)
soft <- pickSoftThreshold(data=expr.t, networkType="signed")
write.table(soft$fitIndices, file="WGCNA_unfilt_thresh.txt",
            col.names=T, row.names=F, quote=F, sep="\t")

library(procona)
unfilt.net <- newProconaObj(networkName="Mtb", 
                     pepdat=expr.t, signed="signed")
#Used power = 4

unfilt.small.net <- newProconaObj(networkName="Mtb", 
                            pepdat=expr.t, signed="signed", minModuleSize=5)

save(unfilt.net, file="WGCNA_net1.RData")
save(unfilt.small.net, file="unfilt.small.net.RData")

#Make a network after filtering for high variance genes
gene.var <- apply(expr.t, 2, function(x) {sd(x, na.rm=T)})


#Plot the variances to look for a good cutoff
plot(sort(gene.var),
     main="Standard Deviation of Gene log2(R/G)",
     xlab="gene",
     ylab="Standard Deviation",
     yaxp=c(0,2.5,10))
grid(col="black")

#Choose 0.75 as cutoff

expr.highvar <- expr.t[,gene.var > 0.75]
high_var.net <- newProconaObj(networkName="Mtb_highvar", 
                            pepdat=expr.highvar, signed="signed")

high_var.thresh <- pickSoftThreshold(expr.highvar)
write.table(high_var.thresh$fitIndices, file="WGCNA_highvar_thresh.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

save(high_var.net, file="high_var.net.RData")

high_var.small <- newProconaObj(networkName="Mtb_highvar", 
                           pepdat=expr.highvar, signed="signed", minModuleSize=5)

#After talking with Shannon, use lower threshold to get more genes
filt_pt5.net <- newProconaObj(networkName="Mtb_pt5_filtered",
                              pepdat=expr.t[,gene.var>0.5], signed="signed",
                              minModuleSize=5)

save(filt_pt5.net, file="filt_pt5.net.RData")

pt5_thresh <- pickSoftThreshold(expr.t[,gene.var > 0.75])

write.table(pt5_thresh$fitIndices, file="WGCNA_pt5_thresh.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

plot(pt5_thresh$fitIndices$Power, pt5_thresh$fitIndices$SFT.R.sq,
     main="2158 Genes SD > 0.5\nSoft R2 vs Power",
     xlab="Power",
     ylab="Soft R2")

plot(pt5_thresh$fitIndices$Power, pt5_thresh$fitIndices$truncated.R.sq,
     main="2158 Genes SD > 0.5\nTruncated R2 vs Power",
     xlab="Power",
     ylab="Soft R2")

write.table(filt_pt5.net@permtest, file="WGCNA_pt5_permtest.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

plot(pickSoftThreshold(filt_pt5.net$fitIndices)


#Plot to see if they follow scale free
scaleFreePlot(unfilt.net@adjacency, main="All genes\n")
scaleFreePlot(high_var.net@adjacency, main="619 Genes SD > 0.75\n")
scaleFreePlot(small_mod@adjacency)
scaleFreePlot(filt_pt5.net@adjacency, main="2158 Genes SD > 0.5\n")
