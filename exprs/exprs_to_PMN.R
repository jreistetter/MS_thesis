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

gse.8786 <- load("GSE8786/g.8786.M.rv.RData")
gpl.4291 <- load("GSE9331/gpl.4291.M.avg.RData")
gpl.4293 <- load("GSE9331/gpl.4293.M.avg.RData")
gpl.8523 <- load("GSE16146/gpl.8523.M.RData")
gpl.8561 <- load("GSE16146/gpl.8561.M.RData")
gpl.8562 <- load("GSE16146/gpl.8562.M.RData")


# 2 - Merge the data frames from each experiment









# 3 - Exclude arrays identified in QC as erroneous
# 4 - Discretize data at a given fold change
# 5 - Write the combined dataframe to a PMN file












