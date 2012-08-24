# Script to look at characteristics of PMN modules for 8.10_hivar including:
#   TF binding
#   Overall module expression

options(stringsAsFactors=F)

#Load data
setwd("~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/")

modules <- read.table("output/8.10_hivar_modules.txt",
                      sep="\t", head=T)
pathways <- read.table("8.10.12_30_mods_0.1_highvar_mods_pathways.txt",
                       sep=":", head=F)