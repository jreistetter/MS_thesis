# Script to take expression data stored as RData objects, discretize them, and 
# then write the data to a file readable by the PMN software.
# 
# 1 - Load expression data from RData objects
# 2 - Discretize platforms without multiple probes
# 3 - Get consensus discretization for platforms with multiple probes
# 4 - Merge the data frames from each experiment
# 5 - Exclude arrays identified in QC as erroneous
# 6 - Filter out arrays with less than 10 DE genes
# 7 - Filter out arrays with more than 10% genes == NA
# 8 - Filter out genes that arent up or down in any arrays
# 9 - Save as txt and RData

# Notes:
#
#  ----- 1.5 fold ----
#   290 initial arrays for express
#  - 24 excluded because of bad MA plots
#  - 23 arrays will be excluded bc < 10 DE genes
#  - 14 arrays excluded bc NA > 10%
#  _____
#   229 total arrays after QC

#   3,842 genes in common between all platforms  
#  -  384 genes that are all no change 
#   _____
#   3,458 genes left after QC

#Fold threshold for discretization
THRESHOLD = 1.5
OUT_FILENAME = "mtb_exprs_1.5fold.txt"
OUT_RFILE = "expr.cons.1.5.RData"

# 1 - Load expression data from RData objects
#

#OHSU wd
#setwd("~/Dropbox")

#Laptop wd
setwd("~/schoolDB/Dropbox")

source("./thesis_work/code/exprs/exprs_funcs.R")

setwd("./thesis_work/data/exprs")

load("GSE8786/g.8786.M.rv.RData")
load("GSE9331/gpl.4291.rv.M.RData")
load("GSE9331/gpl.4293.rv.M.RData")
load("GSE16146/gpl.8523.M.RData")
load("GSE16146/gpl.8561.rv.M.RData")
load("GSE16146/gpl.8562.M.RData")
load("GSE8839/g.8561.rv.M.RData")

#Discretize data sets without multiple probes

g.8786.disc <- discretize(g.8786.M.rv, THRESHOLD)
g.8523.disc <- discretize(gpl.8523.rv.M, THRESHOLD)
g.8562.disc <- discretize(gpl.8562.rv.M, THRESHOLD)

#For datasets with multiple probes, use consensus to roll up to gene


########################################
#
#      GSE8839 GPL 8561 consensus
#
########################################

#Get the total number of unique genes
g.8839.genes <- unique(rownames(g.8561.rv.M))
length(g.8839.genes)
#[1] 3924

#Get the gene IDs of genes with multiple probes
g.8839.dupes <- unique(rownames(g.8561.rv.M)[which(duplicated(rownames(g.8561.rv.M)))])

length(g.8839.dupes)
#[1] 655 genes with multiple probes

#Subset data frame into duplicated and not duplicated
no_dup.8839.idx <- which(!(rownames(g.8561.rv.M) %in% g.8839.dupes))
dup.8839.idx <- which((rownames(g.8561.rv.M) %in% g.8839.dupes))

dup.8839.df <- g.8561.rv.M[dup.8839.idx,]
dim(dup.8839.df)
#[1] 1417  115

#Check that the number of rows of unduplicated genes + rows of duplicated
#genes equals the orignal df size
stopifnot(length(no_dup.8839.idx) + nrow(dup.8839.df) == nrow(g.8561.rv.M))

#Check when we take out duplicated genes we have the number of non-duplicated genes
stopifnot(3924 - 655 == length(no_dup.8839.idx))

cons.8839.df <- df.consensus(dup.8839.df, unique(dup.8839.df$gene), THRESHOLD)
cons.8839.check <- df.check_consensus(dup.8839.df, unique(dup.8839.df$gene), THRESHOLD)
save(cons.8839.check, file="./consensus_check/cons.8839.check.RData")

#Check that it collapsed the genes correctly:
stopifnot(nrow(cons.8839.df) == length(g.8839.dupes))

#discretize the non-duplicated genes
nodup.8839.df <- g.8561.rv.M[no_dup.8839.idx,]
nodup.8839.disc <- discretize(nodup.8839.df, THRESHOLD)
stopifnot(nrow(nodup.8839.disc) == length(no_dup.8839.idx))

#combine the two
g.8839.disc <- rbind(cons.8839.df, nodup.8839.disc[,1:(ncol(nodup.8839.disc))-1])

#Check that the right number of genes are there
stopifnot(nrow(g.8839.disc)==length(g.8839.genes))

#Check that all the genes are there:
stopifnot(length(intersect(rownames(g.8839.disc), g.8839.genes)) == length(g.8839.genes))

#Check that all the arrays are the same
stopifnot(
  sum(colnames(g.8561.rv.M)[1:(ncol(g.8561.rv.M)-1)] != colnames(g.8839.disc)) == 0)

########################################
#
#      GSE9331 GPL 4291 consensus
#
########################################

#Get unique gene IDs
gpl.4291.genes <- unique(rownames(gpl.4291.rv.M))
length(gpl.4291.genes)
#[1] 3897

#Get the gene IDs of genes with multiple probes
gpl.4291.dupes <- unique(rownames(gpl.4291.rv.M)[which(duplicated(rownames(gpl.4291.rv.M)))])
length(gpl.4291.dupes)
#[1] 3897, so all genes have multiple probes

gpl.4291.disc <- df.consensus(gpl.4291.rv.M, gpl.4291.genes, THRESHOLD)
gpl.4291.check <- df.check_consensus(gpl.4291.rv.M, gpl.4291.genes, THRESHOLD)
save(gpl.4291.check, file="./consensus_check/gpl.4291.check.RData")

#Check that all genes are represented
stopifnot(length(intersect(gpl.4291.genes, rownames(gpl.4291.disc)))==length(gpl.4291.genes))

########################################
#
#      GSE9331 GPL 4293 consensus
#
########################################

#Get unique gene IDs
gpl.4293.genes <- unique(rownames(gpl.4293.rv.M))
length(gpl.4293.genes)
#[1] 3900

#Get the gene IDs of genes with multiple probes
gpl.4293.dupes <- unique(rownames(gpl.4293.rv.M)[which(duplicated(rownames(gpl.4293.rv.M)))])
length(gpl.4293.dupes)
#[1] 3900, so all genes have multiple probes

gpl.4293.disc <- df.consensus(gpl.4293.rv.M, gpl.4293.genes, THRESHOLD)
gpl.4293.check <- df.check_consensus(gpl.4293.rv.M, gpl.4293.genes, THRESHOLD)
save(gpl.4293.check, file="./consensus_check/gpl.4293.check.RData")

#Check that all genes are represented
stopifnot(length(intersect(gpl.4293.genes, rownames(gpl.4293.disc)))==length(gpl.4293.genes))


########################################
#
#      GSE16146 GPL 8561 consensus
#
########################################

gpl.8561.rv.M$gene <- rownames(gpl.8561.rv.M)

#Get the total number of unique genes
gpl.8561.genes <- unique(rownames(gpl.8561.rv.M))
length(gpl.8561.genes)
#[1] 3924

#Get the gene IDs of genes with multiple probes
gpl.8561.dupes <- unique(rownames(gpl.8561.rv.M)[which(duplicated(rownames(gpl.8561.rv.M)))])

length(gpl.8561.dupes)
#[1] 655 genes with multiple probes

#Subset data frame into duplicated and not duplicated
no_dup.8561.idx <- which(!(rownames(gpl.8561.rv.M) %in% gpl.8561.dupes))
dup.8561.idx <- which((rownames(gpl.8561.rv.M) %in% gpl.8561.dupes))

dup.8561.df <- gpl.8561.rv.M[dup.8561.idx,]
dim(dup.8561.df)
#[1] 1417  32

#Check that the number of rows of unduplicated genes + rows of duplicated
#genes equals the orignal df size
stopifnot(length(no_dup.8561.idx) + nrow(dup.8561.df) == nrow(gpl.8561.rv.M))

#Check when we take out duplicated genes we have the number of non-duplicated genes
stopifnot(3924 - 655 == length(no_dup.8561.idx))

cons.8561.df <- df.consensus(dup.8561.df, unique(dup.8561.df$gene), THRESHOLD)
cons.8561.check <- df.check_consensus(dup.8561.df, unique(dup.8561.df$gene), THRESHOLD)
save(cons.8561.check, file="./consensus_check/cons.8561.check.RData")

#Check that it collapsed the genes correctly:
stopifnot(nrow(cons.8561.df) == length(gpl.8561.dupes))

#discretize the non-duplicated genes
nodup.8561.df <- gpl.8561.rv.M[no_dup.8561.idx,]
nodup.8561.disc <- discretize(nodup.8561.df, THRESHOLD)
stopifnot(nrow(nodup.8561.disc) == length(no_dup.8561.idx))

#combine the two
gpl.8561.disc <- rbind(cons.8561.df, nodup.8561.disc[,1:(ncol(nodup.8561.disc))-1])

#Check that the right number of genes are there
stopifnot(nrow(gpl.8561.disc)==length(gpl.8561.genes))

#Check that all the genes are there:
stopifnot(length(intersect(rownames(gpl.8561.disc), gpl.8561.genes)) == length(gpl.8561.genes))

#Check that all the arrays are the same
stopifnot(
  sum(colnames(gpl.8561.rv.M)[1:(ncol(gpl.8561.rv.M)-1)] != colnames(gpl.8561.disc)) == 0)


# 2 - Merge the data frames from each experiment

rownames(g.8786.disc) <- toupper(rownames(g.8786.disc))
rownames(gpl.4291.disc) <- toupper(rownames(gpl.4291.disc))
rownames(gpl.4293.disc) <- toupper(rownames(gpl.4293.disc))
rownames(g.8523.disc) <- toupper(rownames(g.8523.disc))
rownames(g.8562.disc) <- toupper(rownames(g.8562.disc))
rownames(gpl.8561.disc) <- toupper(rownames(gpl.8561.disc))
rownames(g.8839.disc) <- toupper(rownames(g.8839.disc))

expr <- merge(g.8786.disc, gpl.4291.disc,
              by.x="row.names", by.y="row.names")

expr <- merge(expr, gpl.4293.disc,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, g.8523.disc,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, g.8562.disc,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, gpl.8561.disc,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, g.8839.disc,
              by.x="Row.names", by.y="row.names")

#Get rid of Row.names column
rownames(expr) <- expr[,1]
expr <- expr[,-1]

#Check that merge worked correctly
stopifnot(
  ncol(expr) == 
    ncol(g.8786.disc) +
    ncol(gpl.4291.disc) +
    ncol(gpl.4293.disc) +
    ncol(g.8523.disc) +
    ncol(g.8562.disc) +
    ncol(gpl.8561.disc) +
    ncol(g.8839.disc)
  )


# 3 - Exclude arrays identified in QC as erroneous

excluded <- read.table("excluded_arrays.txt", head=F, stringsAsFactors=F)[,1]
cols.excl <- which(colnames(expr) %in% excluded)

length(cols.excl)
#[1] 24 excluded because of bad MA plots

#Check that right number of columns will be excluded
stopifnot(length(excluded)==length(cols.excl))


ncol(expr)
#[1] 290
expr <- expr[,-cols.excl]
290 - ncol(expr)
#24, all the arrays excluded

ncol(expr)
#[1] 266

#Filter out arrays with less than 10 DE genes, 9 arrays will be removed
lt.10 <- which(unlist(lapply(expr, function(x) sum(abs(x), na.rm=t))) < 10)
length(lt.10)
#[1] 23 arrays will be excluded bc < 10 DE genes


expr <- expr[,-lt.10]
stopifnot(266 - ncol(expr) == 23)

#Remove any genes that are coded no change for all arrays
all_nochange <- apply(expr, 1, function(x) all(x == 0, na.rm=T))
sum(all_nochange)
#[1] 384 genes that are all no change

#Remove them from the expression dataframe
expr <- expr[!all_nochange,]

dim(expr)
#[1] 3458  243

#Remove any arrays with more than 10% NA values (345 genes)
bad_array <- which(unlist(lapply(expr, function(x) sum(is.na(x)))) > 348)

length(bad_array)
#[1] 14 arrays excluded

expr <- expr[,-bad_array]
dim(expr)
#[1] 3458  229

#check that right number excluded
stopifnot(243-ncol(expr) == 14)

##Write out file
f <- file(OUT_FILENAME, "w")
header <- paste(c(colnames(expr)), collapse="\t")
header <- paste("Name", header, sep="\t")
write(header, f)

write.table(expr, f, quote=F, row.names=T, col.names=F, sep="\t")
close(f)




#save for later use
save(expr, file=OUT_RFILE)
