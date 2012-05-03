# Script for QC and EDA on BUGS expression data

# Written by Joe Reistetter

options(stringsAsFactors=F)

setwd("~/Dropbox/thesis_work/data/")
#Load data
load("exprs/EBUGS58/BUGS58.samples.RData")
load("exprs/EBUGS58/EBUGS58.arrays.RData")

dim(BUGS58.samples)
#[1] 44  4 44 samples

dim(BUGS58.arrays)
#[1] 3765   44  44 samples, 3765 genes

sum(BUGS58.samples$filename != colnames(BUGS58.arrays))
#0, sample df and arrays are in same order

by.cell <- order(BUGS58.samples$celltype)
arr.labels <- BUGS58.samples$celltype[by.cell]
arr.labels[arr.labels=="Aerobic"] <- "A"
arr.labels[arr.labels=="DC"] <- "D"
arr.labels[arr.labels=="MDM"] <- "M"

boxplot(BUGS58.arrays[,by.cell],
        main="E-BUGS-58 Arrays\nnormalized log2(R/G)",
        ylab="normalized log2(R/G)",
        xlab="Arrays",
        names=arr.labels)


perc_DE.15 <- unlist(lapply(BUGS58.arrays, function(x) {
  sum(abs(x) > 1.5, na.rm=T) / sum(!is.na(x))
}))

perc_DE.2 <- unlist(lapply(BUGS58.arrays, function(x) {
  sum(abs(x) > 2, na.rm=T) / sum(!is.na(x))
}))


boxplot(data.frame(fold.1.5=perc_DE.15, fold2=perc_DE.2),
        main="E-BUGS-58 Arrays\nProportion of genes > fold change",
        ylab="Proportion of genes",
        xlab="Fold change threshold")

boxplot(data.frame(fold.1.5=perc_DE.15, fold2=perc_DE.2),
        main="E-BUGS-58 Arrays\nProportion of genes > fold change",
        ylab="Proportion of genes",
        xlab="Fold change threshold")


boxplot(perc_DE.2 ~ as.factor(BUGS58.samples$celltype),
        main="E-BUGS-58 Arrays\nProportion of genes > 2 fold change by cell type",
        ylab="Proportion of genes",
        xlab="Cell type")


