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
# Basic algorithm to assign p-DNA edges:
#   1. Make a dataframe, each row is an Rv ID and its operon
#   2. Make a list, with element IDs corresponding to operon IDs, 
#     and the element value is a character vector of Rv IDs
#   3. For each regulator, look up its targets' operon IDs
#   4. Look up operon members in list
#   5. Assign a p-DNA edge between the regulator and all the operon members.
#
# Results:
#   -DOOR has 917 operons with 2607 genes as door_ops, by genes in door_genes
#   -Microbes Online has 496 operons with 1259 genes as microbes_ops, 
#                     by genes in microbes_genes
#
#
# Written by Joe Reistetter

# -import data

#At OHSU Dropbox
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
setwd("~/schoolDB/Dropbox/thesis_work/data/protein-DNA/")

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
door_ops <- list()

for (operon in unique(door_genes$OperonID)){
  op_genes <- door_genes[door_genes$OperonID == operon,2]
  door_ops[[as.character(operon)]] <- op_genes
}


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

#Initializes variables to use in loop that creates operons
microbes_ops <- list()
i <- 1
op_id <- 1

#Loop through each row of the filtered operon predictions.
#Depending on the continuity of the operon, loop does
#different things.
while (i < nrow(mic_p_filtered)+1){
  
  #If there is an NA, skip that row
  if(is.na(mic_p_filtered[i+1,1]) | is.na(mic_p_filtered[i,2])){
    i <- i+1
    op_id <- op_id+1
    next
  }
    
  #If the operon is only 2 genes, add those 2 genes as a new operon
  if(mic_p_filtered[i+1,1] != mic_p_filtered[i,2]){
    op_id <- op_id+1
    microbes_ops[[as.character(op_id)]] <- c(microbes_ops[[as.character(op_id)]],
                                             unlist(c(mic_p_filtered[i,c(3,4)])))
    i <- i+1
    op_id <- op_id+1
    next
  }
  
  #While the genes are continuous (in an operon), keep addding to new operon
  while(mic_p_filtered[i+1,1] == mic_p_filtered[i,2]){

    existing <- microbes_ops[[as.character(op_id)]]
    these_genes <- unlist(c(mic_p_filtered[i,c(3,4)], mic_p_filtered[i+1,c(3,4)]))
    genes_add <- c(existing, these_genes)
    microbes_ops[[as.character(op_id)]] <- genes_add
    i <- i+1
  }
  op_id <- op_id+1
  i <- i+1
}

#Add last operon
microbes_ops[[length(microbes_ops)+1]] <- 
  unlist(c(mic_p_filtered[i-1,c(3,4)]))

#For simplicity, the while loop doesn't check if a gene is already in an operon
#before adding it. This removes duplicates. Based on algorithm, duplicates
#are expected..,.,.. 

microbes_ops <- lapply(microbes_ops, unique)

microbes_genes <- data.frame(operon_ID=c(), gene=c())

for (i in c(1:length(microbes_ops))){
  genes <- microbes_ops[[i]]
  rows <- matrix(c(rep(i, length(genes)), genes), ncol=2)
  microbes_genes <- rbind(microbes_genes, rows)
}

colnames(microbes_genes) <- c("operon_ID", "gene_ID")


#Check to see if all the filtered gene IDs are in the operons
sum(!(c(mic_p_filtered$SysName1, mic_p_filtered$SysName2) 
      %in% microbes_genes$gene_ID))
#0, yes

###################################################
#
#           TBDB data
#
###################################################

#Not sure if going to use this operon data yet since not sure how it was made.
#Keep code around for assigning regulators to targets.


#OHSU WD
#setwd("~/Dropbox/thesis_work/data/ChIP-Seq/")
#laptop WD
setwd("~/schoolDB/Dropbox/thesis_work/data/ChIP-Seq/")


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