# Script to take expression data stored as RData objects, discretize them, and 
# then write the data to a file readable by the PMN software.
# 
# 1 - Load expression data from RData objects
# 2 - Merge the data frames from each experiment
# 3 - Exclude arrays identified in QC as erroneous
# 4 - Discretize data at a given fold change
# 5 - Write the combined dataframe to a PMN file

# 1 - Load expression data from RData objects
#

#OHSU wd
#setwd("~/Dropbox)

#Laptop wd
setwd("~/schoolDB/Dropbox")

setwd("./thesis_work/data/exprs")

load("GSE8786/g.8786.M.rv.RData")
load("GSE9331/gpl.4291.M.avg.RData")
load("GSE9331/gpl.4293.M.avg.RData")
load("GSE16146/gpl.8523.M.RData")
#Dont use 8561 until duplication issue is resolved
#load("GSE16146/gpl.8561.M.RData")
load("GSE16146/gpl.8562.M.RData")


# 2 - Merge the data frames from each experiment

rownames(g.8786.M.rv) <- toupper(rownames(g.8786.M.rv))
rownames(gpl.4291.M.avg) <- toupper(rownames(gpl.4291.M.avg))
rownames(gpl.4293.M.avg) <- toupper(rownames(gpl.4293.M.avg))
rownames(gpl.8523.rv.M) <- toupper(rownames(gpl.8523.rv.M))
rownames(gpl.8562.rv.M) <- toupper(rownames(gpl.8562.rv.M))

expr <- merge(g.8786.M.rv, gpl.4291.M.avg,
              by.x="row.names", by.y="row.names")

expr <- merge(expr, gpl.4293.M.avg,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, gpl.8523.rv.M,
              by.x="Row.names", by.y="row.names")

expr <- merge(expr, gpl.8562.rv.M,
              by.x="Row.names", by.y="row.names")

#save for later use
save(expr, "expr.RData")

#get rid of all objects except expr
objs <- ls()
objs[which(objs == "expr")] <- "objs"

rm(list=objs)

# 3 - Exclude arrays identified in QC as erroneous
# 4 - Discretize data at a given fold change
# 5 - Write the combined dataframe to a PMN file














expr.objs <- ls()

mmerge <- function(obj.names){
  command <- sprintf('merge(%s, %s, by.x="row.names", by.y="row.names")', obj.names[1], obj.names[2])
  df <- eval(parse(text=command))
  for (i in c(3:length(obj.names))){
    command <- sprintf('merge(df, %s, by.x="row.names", by.y="row.names")', obj.names[i])
    df <- eval(parse(text=command))
  }
  
  return(df)
}

expr <- mmerge(expr.objs)


