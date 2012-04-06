# Script to analyze the coefficient of variation for probes
# with technical replicates
# 
# Written by Joe Reistetter

#ohsu
db <- "~/Dropbox/thesis_work/"
#laptop
#db <- "~/schoolDB/Dropbox/thesis_work/"

source(paste(db, "code/exprs/exprs_funcs.R", sep=""))

######################################################
#
#         GPL 4291
#
#
######################################################

path.4291 <- paste(db, "data/exprs/GSE9331/gpl.4291.bc.norm.rv.RData", sep="")

load(path.4291)

gpl.4291.rv.M <- remove_bad_spots(gpl.4291.bc.norm.rv)
gpl.4291.rv.A <- remove_bad_spots(gpl.4291.bc.norm.rv, type="A")


#Get coefficient of variation of ratios of the the replicated probes
#Remove bad arrays first
bad.4291 <- c("GSM237637.gpr", "GSM237644.gpr", "GSM237645.gpr", "GSM237646.gpr", "GSM237650.gpr", 
              "GSM237651.gpr", "GSM237655.gpr", "GSM237666.gpr", "GSM237667.gpr")

good.4291.col.idx <- which(!(colnames(gpl.4291.rv.M) %in% bad.4291))
gpl.4291.M <- gpl.4291.rv.M[,good.4291.col.idx]
gpl.4291.A <- gpl.4291.rv.A[,good.4291.col.idx]

gpl.4291.M$gene <- rownames(gpl.4291.M)
gpl.4291.A$gene <- rownames(gpl.4291.A)

replicated.4291 <- unique(gpl.4291.M[duplicated(gpl.4291.M$gene),]$gene)
coefs.4291 <- sapply(replicated.4291, probe_CV, df=gpl.4291.M)
mean.A.4291 <- sapply(replicated.4291, probe_mean_A, df=gpl.4291.A)
mean.M.4291 <- sapply(replicated.4291, probe_mean_A, df=gpl.4291.M)


#########
#
#     Plots of CV vs M and A
#
#########

#CV probably isn't a good measure to filter genes with technical replicate probes,
#these plots will show that
par(mfrow=c(2,2))

#Plots of all the CVs vs M and A
plot(unlist(mean.A.4291), unlist(coefs.4291),
     main='CV as a function of A\nGPL4291',
     xlab='1/2 * log2(R*G)',
     ylab='CV')

#This plot shows that the extreme values are all near 0
#which makes sense given CV = SD/mean
plot(unlist(mean.M.4291), unlist(coefs.4291),
     main = 'CV as a function of log2 (R/G)\nGPL4291',
     xlab='log2(R/G)',
     ylab = 'CV')

#Filter out extreme values to see if this trend is true across the data
no_outliers.4291 <- which(abs(unlist(coefs.4291))<10)

plot(unlist(mean.A.4291)[no_outliers.4291], unlist(coefs.4291)[no_outliers.4291],
     main = 'CV as a function of A\nFiltered CV < 10\nGPL4291',
     xlab='1/2 * log2(R*G)',
     ylab = 'CV')

plot(unlist(mean.M.4291)[no_outliers.4291], unlist(coefs.4291)[no_outliers.4291],
     main = 'CV as a function of log2 (R/G)\nFiltered CV < 10\nGPL4291',
     xlab='log2(R/G)',
     ylab = 'CV')

coefs.mean.4291 <- apply(coefs.4291, 2, mean)
coefs.median.4291 <- apply(coefs.4291, 2, median)

#Count the number of genes that would be excluded at CoV < 1
excl.4291 <- sum(unlist(lapply(coefs.4291, function(x) sum(x > 1, na.rm=T)))) #20,388 excluded

excl.4291.2 <- sum(unlist(lapply(coefs.4291, function(x) sum(x > 2, na.rm=T)))) #12,147 excluded

#Total number of genes:
dim(coefs.4291)[1] * dim(coefs.4291)[2] 
#120,807, so roughly 20% would be excluded at 1 threshold

#Look at CoV on a per-gene basis to see if some gene are all bad

gene.cov.4291 <- apply(coefs.4291, 2, function(x) sum(x > 1, na.rm=T))
max(gene.cov.4291) #17


png("GPL4291_probe_mean_CoV.png")
boxplot(coefs.mean.4291, 
        main = "Gene-wise mean coefficient of variation\nProbe log(R/G)\nGPL4291 22 arrays 3897 genes",
        ylab="Mean CoV")
dev.off()

png("GPL4291_probe_median_CoV.png")
boxplot(coefs.median.4291, 
        main = "Gene-wise median coefficient of variation\nProbe log(R/G)\nGPL4291 22 arrays 3897 genes",
        ylab="Median CoV")
dev.off()








######################################################
#
#         GPL 4293
#
#
######################################################

replicated.4293 <- unique(gpl.4293.rv.M[duplicated(gpl.4293.rv.M$gene),]$gene)
coefs.4293 <- sapply(replicated.4293, probe_CV, df=gpl.4293.rv.M)
coefs.mean.4293 <- apply(coefs.4293, 2, mean)
coefs.median.4293 <- apply(coefs.4293, 2, median)

#Count the number of genes that would be excluded at CoV < 1
excl.4293 <- sum(unlist(lapply(coefs.4293, function(x) sum(x > 1, na.rm=T)))) #3595

excl.4293.2 <- sum(unlist(lapply(coefs.4293, function(x) sum(x > 2, na.rm=T)))) #1823 excluded

#Total number of genes:
dim(coefs.4293)[1] * dim(coefs.4293)[2] 
#42,900, so roughly 10% would be excluded at 1 threshold

#Look at CoV on a per-gene basis to see if some gene are all bad

gene.cov.4293 <- apply(coefs.4293, 2, function(x) sum(x > 1, na.rm=T))
max(gene.cov.4293) #9

png("GPL4293_probe_mean_CoV.png")
boxplot(coefs.mean.4293, 
        main = "Gene-wise mean coefficient of variation\nProbe log(R/G)\nGPL4293 11 arrays 3897 genes",
        ylab="Mean CoV")
dev.off()

png("GPL4293_probe_median_CoV.png")
boxplot(coefs.median.4293, 
        main = "Gene-wise median coefficient of variation\nProbe log(R/G)\nGPL4293 11 arrays 3897 genes",
        ylab="Median CoV")
dev.off()