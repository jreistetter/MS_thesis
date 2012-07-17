# Script to create input file to PMN of only highly variable genes, using
# the same set that was used to create the WGCNA network.
# 
# Written by Joe Reistetter

setwd("~/Dropbox/thesis_work/data/exprs/")
load("expr.cons.1.5.RData")
load("filt_pt5.net.RData")

#Grab the genes used in building the WGCNA network
wgcna_genes <- filt_pt5.net@peptides

length(wgcna_genes)
# [1] 2158 genes

#Filter the discretized expression to WGCNA genes
expr.wgcna <- expr[which(rownames(expr) %in% wgcna_genes), ]

#Write out in PMN format
f <- file("mtb_high_var_expr_1.5.txt", "w")
header <- paste(c(colnames(expr.wgcna)), collapse="\t")
header <- paste("Name", header, sep="\t")
write(header, f)

write.table(expr.wgcna, f, quote=F, row.names=T, col.names=F, sep="\t")
close(f)