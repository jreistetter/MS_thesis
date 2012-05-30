# Script to visualize the WGCNA modules and coexpression network
# 
# Written by Joe Reistetter

library(WGCNA)

#load data
setwd("~/Dropbox/thesis_work/")

load("data/exprs/filt_pt5.net.RData")

plotDendroAndColors(filt_pt5.net@geneTree, filt_pt5.net@mergedColors,
                    "Modules",
                    dendroLabels=F,
                    hang=0.03,
                    addGuide=T,
                    guideHang=0.05,
                    main="")

