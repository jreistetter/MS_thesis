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
#   -DOOR has 917 operons
# Written by Joe Reistetter

# -import data

#At OHSU Dropbox
#setwd("/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/")
#My laptop Dropbox
setwd("~/schoolDB/Dropbox/thesis_work/data/protein-DNA/")

###################################################
#
#DOOR data
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
#Microbes Online data
#
###################################################

#Note: pOp column is the probability that the two are in an operon together
#     0 = not at all, 1 = absolutely in the same operon

mic_online <- read.table("microbes_online_mtb_operons.txt", head=T,
                         stringsAsFactor=F, sep="\t", quote="")
mic_p_filtered <- mic_online[which(mic_online$pOp >= 0.9),]

microbes_ops <- list()
i <- 1
op_id <- 1

while (i < nrow(mic_p_filtered)+1){
  if(is.na(mic_p_filtered[i+1,1]) | is.na(mic_p_filtered[i,2])){
    i <- i+1
    op_id <- op_id+1
    next
  }
    
  
  if(mic_p_filtered[i+1,1] != mic_p_filtered[i,2]){
    op_id <- op_id+1
    microbes_ops[[as.character(op_id)]] <- c(microbes_ops[[as.character(op_id)]],
                                             unlist(c(mic_p_filtered[i,c(3,4)])))
    i <- i+1
    op_id <- op_id+1
    next
  }
  
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

microbes_ops <- lapply(microbes_ops, unique)


# -load H37Rv annotation obtained from tbdb.org
h37rv.annot <- read.table("tbdb_H37rv_annotation.txt", sep="\t", head=T, stringsAsFactors=F,
                          colClasses=c(rep("NULL",3), "numeric", "numeric", "character",
                                       rep("NULL",3), "character", rep("NULL",3)),
                          quote=""
                          )

nrow(h37rv.annot) == 3999 #3999 is how many genes tbdb.org said were downloaded

setwd("../protein-DNA/")

# -use gene start/stop and operon start/stop coords to assign genes to operons
#  do minus and plus strands separate
library(IRanges)

# Minus strand
operons.minus <- operons[operons$Strand=="-",]
genes.minus <- h37rv.annot[h37rv.annot$Strand=="-",]
rownames(genes.minus) <- c(1:nrow(genes.minus))

o.minus.r <- RangedData(ranges=IRanges(start = operons.minus$Start,
                                    end = operons.minus$Stop),
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
operons.plus <- operons[operons$Strand=="+",]
genes.plus <- h37rv.annot[h37rv.annot$Strand=="+",]
rownames(genes.plus) <- c(1:nrow(genes.plus))

o.plus.r <- RangedData(ranges=IRanges(start = operons.plus$Start,
                                      end = operons.plus$Stop),
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
nrow(operons) == (nrow(operons.minus) + nrow(operons.plus))


#Now need to assign regulators to operons


#  making sure that the gene and operon are on the same strand
# -for each gene in an operon, create an interaction between regulator and gene (protein-DNA edge)