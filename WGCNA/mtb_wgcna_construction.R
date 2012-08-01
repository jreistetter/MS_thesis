# Discover modules via WGCNA using a power of 8 for the soft-thresholding
# in the adjacency matrix and save as an R object for use in later scripts.
# 
# Written by Joe Reistetter
# Last updated: 7/31/2012
#
# -Power is 0.6
# -Minimum module size is 10
# -Genes filtered for SD > 0.5
# -TOM is signed

library(WGCNA)
library(procona)

# Load the expression data (fold change)
setwd("~/Dropbox/thesis_work/")
load("./data/exprs/expr.mean.RData")

# Calculate standard deviation to filter genes. See expr_wgcna.R for
# plots that show how cutoff of 0.5 SD was reached.
gene.var <- apply(expr.mean, 1, function(x) {sd(x, na.rm=T)})

# Filter the genes to retain ones with high variation
expr.highvar <- expr.mean[gene.var > 0.5,]

# Transpose to match WGCNA input format
expr.highvar.t <- t(expr.highvar)

# Get stats to pick soft threshold
thresholds <- pickSoftThreshold(expr.highvar.t)

# Choose power of 6 (r^2 = 0.88, trunc=0.987):
# mean.k. median.k. max.k.
# 7.9900  4.44      54.00

mtb.procona.6 <- newProconaObj("mtb_6", #network name
                          pepdat=expr.highvar.t, #expression data
                          pow=6, #set power to 6 for soft threshold
                          signed="signed", #sign considered for adjacency
                          minModuleSize=10)

save(mtb.procona.6, file="data/WGCNA/mtb.procona.6.RData")

# Visualize module dendrogram and TOM pairwise heatmap

plotDendroAndColors(mtb.procona.6@geneTree, mtb.procona.6@mergedColors,
                    "Modules",
                    dendroLabels=F,
                    hang=0.03,
                    addGuide=T,
                    guideHang=0.05,
                    main="")


dissTOM <- 1 - mtb.procona.6@TOM
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
TOMplot(plotTOM, mtb.procona.6@geneTree, mtb.procona.6@mergedColors)