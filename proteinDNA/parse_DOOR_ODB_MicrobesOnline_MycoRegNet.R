# Script to take operon raw data from DOOR, Microbes Online, ODB, MycoRegNet,
# and parse into R objects for use in assigning protein-DNA edges for PMN.
# 
# Also parses in regulator-target pairs from MycoRegNet (article)
#
# Results:
#   -DOOR has 917 operons with 2607 genes as door_ops, by genes in door_genes
#   -Microbes Online has 496 operons with 1259 genes as microbes_ops, 
#                     by genes in microbes_genes
#   -ODB had 39 operons with 99 genes as ODB.op.list, ODB_genes
#   -MycoRegNet has 180 reg-target pairs from TF predictions, 223 from
#   orthologous predictions from other species.
#   
# Written by Joe Reistetter

# -import data

#At OHSU Dropbox
root_path <-"~/Dropbox/thesis_work"

setwd(paste(root_path, "/data/protein-DNA/", sep=""))
source(paste(root_path, "/code/thesis_funcs.R", sep=""))

###################################################
#
#                 DOOR data
#
###################################################

#Read in data
door <- read.table("door_operons.txt", head=T, stringsAsFactors=F,
                   sep="\t", quote="")

#Extract gene-operon assignments
door_genes <- door[,c(1,3)]

#Store operon members in a list
door.op.list <- list()

for (operon in unique(door_genes$OperonID)){
  op_genes <- door_genes[door_genes$OperonID == operon,2]
  door.op.list[[as.character(operon)]] <- op_genes
}

save(door.op.list, file="./RData/door.op.list.RData")
save(door_genes, file="./RData/door_genes.RData")

###################################################
#
#           Microbes Online data
#
###################################################

#Note: pOp column is the probability that the two are in an operon together
#     0 = not at all, 1 = absolutely in the same operon

#Load data
mic_online <- read.table("microbes_online_mtb_operons.txt", head=T,
                         stringsAsFactor=F, sep="\t", quote="")

#Filter out operon predictions with pOp >= 0.9
mic_p_filtered <- mic_online[which(mic_online$pOp >= 0.9),]

microbes.op.list <- operons_from_gene_pairs(mic_p_filtered, c(3,4))

microbes_genes <- op_list_to_df(microbes.op.list)

#Check to see if all the filtered gene IDs are in the operons
sum(!(c(mic_p_filtered$SysName1, mic_p_filtered$SysName2) 
      %in% microbes_genes$gene_ID))
#0, yes

save(microbes.op.list, file="./RData/microbes.op.list.RData")
save(microbes_genes, file="./RData/microbes_genes.RData")

###################################################
#
#           ODB Data
#
###################################################

ODB <- read.table("ODB_operons.txt", head=T, stringsAsFactors=F,
                  sep="\t", quote="")

#genes column contains a space delimited string of operon members
#parse into list
ODB.op.list <-  list()

for (i in 1:nrow(ODB)){
  genes.str <- ODB[i,]$genes
  genes <- unlist(strsplit(genes.str, " ", fixed=T))
  ODB.op.list[[as.character(i)]] <- genes
}

ODB_genes <- op_list_to_df(ODB.op.list)

save(ODB.op.list, file="./RData/ODB.op.list.RData")
save(ODB_genes, file="./RData/ODB_genes.RData")

###################################################
#
#           MycoRegNet Data (regulators and targets)
#
###################################################

#Table 1 predicts regulators from orthologous genes in other species
myco.pairs.ortho <- read.table("./MycoRegNet/mycoregnet_tbl1_parsed.txt",
                         head=F, sep="\t", stringsAsFactors=F,
                         col.names=c("regulator", "target"))

#Table 2 predicts regulators from conserved binding sites
myco.pairs.TFBS<- read.table("./MycoRegNet/mycoregnet_tbl2_parsed.txt",
                         head=F, sep="\t", stringsAsFactors=F,
                         col.names=c("regulator", "target"))


save(myco.pairs.ortho, file="./RData/myco.pairs.ortho.RData")
save(myco.pairs.TFBS, file="./RData/myco.pairs.TFBS.RData")


###################################################
#
#           TBDB data (not used but saving code)
#
###################################################

#Not sure if going to use this operon data yet since not sure how it was made.
#Keep code around for assigning regulators to targets.

setwd(paste(root_path, "/data/ChIP-Seq/", sep=""))


tbdb.reg_targ <- read.table("peaks_2012-01-06_18-35-16.tsv", sep="\t",
                       head=T, stringsAsFactors=F, 
                       colClasses=c("character", "character", "NULL","NULL","NULL",
                                    "numeric", "numeric", "NULL","NULL","NULL")
                       )
head(tbdb.reg_targ)

tbdb.operons <- read.table("tbdb_operons.txt", sep=",", head=T, stringsAsFactors=F,
                      colClasses=c(rep("NULL",3), "numeric", "numeric",
                                   "character", rep("NULL",5)),
                      quote=""
                      )

nrow(tbdb.operons) == 2525 #Number of operons downloaded from tbdb.org

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
tbdb.operons.minus <- tbdb.operons[tbdb.operons$Strand=="-",]
genes.minus <- h37rv.annot[h37rv.annot$Strand=="-",]
rownames(genes.minus) <- c(1:nrow(genes.minus))

o.minus.r <- RangedData(ranges=IRanges(start = tbdb.operons.minus$Start,
                                    end = tbdb.operons.minus$Stop),
                                    space = 1)
g.minus.r <- RangedData(ranges=IRanges(start = genes.minus$Start,
                                       end = genes.minus$Stop),
                        space = 1)

minus.overlaps <- as.data.frame(as.matrix(findOverlaps(g.minus.r, o.minus.r)))
colnames(minus.overlaps) <- c("gene.idx", "operon")



#add the locus ID
minus.overlaps <- merge(minus.overlaps, genes.minus, by.x="gene.idx", by.y="row.names",
                        all.x=T)
minus.overlaps <- minus.overlaps[with(minus.overlaps, order(operon)),]


# Positive strand
tbdb.operons.plus <- tbdb.operons[tbdb.operons$Strand=="+",]
genes.plus <- h37rv.annot[h37rv.annot$Strand=="+",]
rownames(genes.plus) <- c(1:nrow(genes.plus))

o.plus.r <- RangedData(ranges=IRanges(start = tbdb.operons.plus$Start,
                                      end = tbdb.operons.plus$Stop),
                       space = 1)
g.plus.r <- RangedData(ranges=IRanges(start = genes.plus$Start,
                                      end = genes.plus$Stop),
                       space = 1)

plus.overlaps <- as.data.frame(as.matrix(findOverlaps(g.plus.r, o.plus.r)))
colnames(plus.overlaps) <- c("gene.idx", "operon")
plus.overlaps <- plus.overlaps[with(plus.overlaps, order(operon)),]

#add the locus ID
plus.overlaps <- merge(plus.overlaps, genes.plus, by.x="gene.idx", by.y="row.names",
                       all.x=T)
plus.overlaps <- plus.overlaps[with(plus.overlaps, order(operon)),]

get_op_regulators <- function(operon, df, reg_df){
  genes <- df[df$operon==operon,6]
  regulators <- c()
  for (i in 1:length(genes)){
    gene_regs <- reg_df[reg_df$target == genes[i], 1]
    regulators <- c(regulators, gene_regs)
  }
  return(unique(regulators))
}



#Make sure that the splits worked
nrow(h37rv.annot) == (nrow(genes.minus) + nrow(genes.plus))
nrow(tbdb.operons) == (nrow(tbdb.operons.minus) + nrow(tbdb.operons.plus))


#Now need to assign regulators to operons


#  making sure that the gene and operon are on the same strand
# -for each gene in an operon, create an interaction between regulator and gene (protein-DNA edge)