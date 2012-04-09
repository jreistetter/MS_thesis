# Script to take expression data stored as RData objects, discretize them, and 
# then write the data to a file readable by the PMN software.
# 
# 1 - Load expression data from RData objects
# 2 - Merge the data frames from each experiment
# 3 - Exclude arrays identified in QC as erroneous
# 4 - Discretize data at a fold thresholds {1, 1.5, 2}
#     -Filter out arrays with less than 10 DE genes
# 5 - Write the combined dataframe to a PMN file

# Notes:
#   128 total arrays for express
#   119 at 1-fold threshold
#   105 at 1.5-fold threshold
#   97 at 2-fold threshold
#   *After filtering for arrays with < 10 DE genes



# 1 - Load expression data from RData objects
#

#OHSU wd
setwd("~/Dropbox")

#Laptop wd
#setwd("~/schoolDB/Dropbox")

setwd("./thesis_work/data/exprs")

load("GSE8786/g.8786.M.rv.RData")
load("GSE9331/gpl.4291.M.avg.RData")
load("GSE9331/gpl.4293.M.avg.RData")
load("GSE16146/gpl.8523.M.RData")
#Dont use 8561 until duplication issue is resolved
#load("GSE16146/gpl.8561.M.RData")
load("GSE16146/gpl.8562.M.RData")
load("GSE8839/g.8561.rv.M.avg.RData")


# 2 - Merge the data frames from each experiment

rownames(g.8786.M.rv) <- toupper(rownames(g.8786.M.rv))
rownames(gpl.4291.M.avg) <- toupper(rownames(gpl.4291.M.avg))
rownames(gpl.4293.M.avg) <- toupper(rownames(gpl.4293.M.avg))
rownames(gpl.8523.rv.M) <- toupper(rownames(gpl.8523.rv.M))
rownames(gpl.8562.rv.M) <- toupper(rownames(gpl.8562.rv.M))
rownames(g.8561.rv.M.avg) <- toupper(rownames(g.8561.rv.M.avg))

expr <- merge(g.8786.M.rv, gpl.4291.M.avg,
              by.x="row.names", by.y="row.names")

expr <- merge(expr, gpl.4293.M.avg,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, gpl.8523.rv.M,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, gpl.8562.rv.M,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, g.8561.rv.M.avg,
              by.x="Row.names", by.y="row.names")

#Get rid of Row.names column
rownames(expr) <- expr[,1]
expr <- expr[,-1]

ncol(expr) == ncol(g.8786.M.rv) + ncol(g.8561.rv.M.avg) + ncol(gpl.4291.M.avg) + ncol(gpl.4293.M.avg) + ncol(gpl.8523.rv.M) + ncol(gpl.8562.rv.M)

#get rid of all objects except expr
objs <- ls()
objs[which(objs == "expr")] <- "objs"

rm(list=objs)

# 3 - Exclude arrays identified in QC as erroneous

excluded <- read.table("excluded_arrays.txt", head=F, stringsAsFactors=F)[,1]
cols.excl <- which(colnames(expr) %in% excluded)
length(excluded)
length(cols.excl)
#List of excluded arrays is 24, since we aren't doing 
#gpl.8561 right now, 2 of the excluded arrays aren't in expr. So discrepancy
#ok.

ncol(expr)
#258
expr <- expr[,-cols.excl]
258 - ncol(expr)
#22, all the arrays excluded


#save for later use
save(expr, file="expr.RData")


# 4 - Discretize data at a given fold change
source("../../code/exprs/exprs_funcs.R")

expr.1 <- discretize(expr, 1)
ncol(expr.1)
#128 arrays

#Filter out arrays with less than 10 DE genes, 9 arrays will be removed
lt.10 <- which(unlist(lapply(expr.1, function(x) sum(abs(x), na.rm=t))) < 10)
expr.1 <- expr.1[,-lt.10]
128 - ncol(expr.1)
#9 arrays removed

expr.1.5 <- discretize(expr, 1.5)
ncol(expr.1.5)
#128 arrays

#Filter out arrays with less than 10 DE genes, 9 arrays will be removed
lt.10 <- which(unlist(lapply(expr.1.5, function(x) sum(abs(x), na.rm=t))) < 10)
expr.1.5 <- expr.1.5[,-lt.10]
128 - ncol(expr.1.5)
#23 arrays removed

expr.2 <- discretize(expr, 2)
ncol(expr.2)
#128 arrays

#Filter out arrays with less than 10 DE genes, 9 arrays will be removed
lt.10 <- which(unlist(lapply(expr.2, function(x) sum(abs(x), na.rm=t))) < 10)
expr.2 <- expr.2[,-lt.10]
128 - ncol(expr.2)
#31 arrays removed

# 5 - Write the combined dataframe to a PMN file


##Write out 1-fold file
f.1 <- file("mtb_exprs_1fold.txt", "w")
header <- paste(c(colnames(expr.1)), collapse="\t")
header <- paste("Name", header, sep="\t")
write(header, f.1)

write.table(expr.1, f.1, quote=F, row.names=T, col.names=F, sep="\t")
close(f.1)


##Write out 1.5 fold file
f.1.5 <- file("mtb_exprs_1.5fold.txt", "w")
header <- paste(c(colnames(expr.1.5)), collapse="\t")
header <- paste("Name", header, sep="\t")
write(header, f.1.5)

write.table(expr.1.5, f.1.5, quote=F, row.names=T, col.names=F, sep="\t")
close(f.1.5)


##Write out 2 fold file
f.2 <- file("mtb_exprs_2fold.txt", "w")
header <- paste(c(colnames(expr.2)), collapse="\t")
header <- paste("Name", header, sep="\t")
write(header, f.2)

write.table(expr.2, f.2, quote=F, row.names=T, col.names=F, sep="\t")
close(f.2)


