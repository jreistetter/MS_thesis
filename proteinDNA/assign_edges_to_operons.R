# Script to take regulator-target lists and operon data from tbdb.org
# and parse it into a list of protein-DNA interactions for use in the PMN software.
# 
# The regulator-target data has the start and stop coords of the peak
# 
# The operon data has:
#   the start and stop coords of the operon
#   the strand of the operon
#   the operon name, which isn't very useful for determining membership
#   the length
# 
# Written by Joe Reistetter

# -import both data

#At OHSU Dropbox
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
setwd("~/schoolDB/Dropbox/thesis_work/data/ChIP-Seq/")
reg_targ <- read.table("peaks_2012-01-06_18-35-16.tsv", sep="\t",
                        head=T, stringsAsFactors=F, 
                       colClasses=c("character", "character", "NULL","NULL","NULL",
                                      "numeric", "numeric", "NULL","NULL","NULL")
                       )
head(reg_targ)

operons <- read.table("tbdb_operons.txt", sep=",", head=T, stringsAsFactors=F,
                      colClasses=c(rep("NULL",3), "numeric", "numeric",
                                   "character", rep("NULL",5)),
                      quote=""
                      )

nrow(operons) == 2525 #Number of operons downloaded from tbdb.org

# -load H37Rv annotation obtained from tbdb.org
h37rv.annot <- read.table("tbdb_H37rv_annotation.txt", sep="\t", head=T, stringsAsFactors=F,
                          colClasses=c(rep("NULL",3), "numeric", "numeric", "character",
                                       rep("NULL",3), "character", rep("NULL",3)),
                          quote=""
                          )

nrow(h37rv.annot) == 3999 #3999 is how many genes tbdb.org said were downloaded

# -use gene start/stop and operon start/stop coords to assign genes to operons
#  do minus and plus strands separate
library(IRanges)

# Minus strand
operons.minus <- operons[operons$Strand=="-",]
genes.minus <- h37rv.annot[h37rv.annot$Strand=="-",]

o.minus.r <- RangedData(ranges=IRanges(start = operons.minus$Start,
                                    end = operons.minus$Stop),
                                    space = 1)
g.minus.r <- RangedData(ranges=IRanges(start = genes.minus$Start,
                                       end = genes.minus$Stop),
                        space = 1)

minus.overlaps <- as.data.frame(as.matrix(findOverlaps(g.minus.r, o.minus.r)))
colnames(minus.overlaps) <- c("gene", "operon")
minus.overlaps <- minus.overlaps[with(minus.overlaps, order(operon)),]

# Positive strand
operons.plus <- operons[operons$Strand=="+",]
genes.plus <- h37rv.annot[h37rv.annot$Strand=="+",]

o.plus.r <- RangedData(ranges=IRanges(start = operons.plus$Start,
                                      end = operons.plus$Stop),
                       space = 1)
g.plus.r <- RangedData(ranges=IRanges(start = genes.plus$Start,
                                      end = genes.plus$Stop),
                       space = 1)

plus.overlaps <- as.data.frame(as.matrix(findOverlaps(g.plus.r, o.plus.r)))
colnames(plus.overlaps) <- c("gene", "operon")
plus.overlaps <- plus.overlaps[with(plus.overlaps, order(operon)),]

#Make sure that the splits worked
nrow(h37rv.annot) == (nrow(genes.minus) + nrow(genes.plus))
nrow(operons) == (nrow(operons.minus) + nrow(operons.plus))


#Now need to assign regulators to operons

dd[with(dd, order(-z, b)), ]

#  making sure that the gene and operon are on the same strand
# -for each gene in an operon, create an interaction between regulator and gene (protein-DNA edge)