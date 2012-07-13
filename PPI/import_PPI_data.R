# Script to create the protein-protein edge list for input to learn a PMN. 
# 
# 1. Import protein-protein interactions from two datasets:
#   -STRING Mtb H37Rv functional protein associations from flatfile
#   -Protein-protein interactions derived from bacterial 2-hybrid assays by Wang et al. 2010 PMID: 20973567
# 
# 2. Filter STRING associations for confidence > 0.9
#
# () May need to use the second supp table from paper to filter
# 
# 3. Calculate the confidence score for each edge as given in Novershtern et al. 2011 PMID: 21685068 
# 
# 4. Write the edges and their confidence scores to a file readable by the PMN software from http://www.compbio.cs.huji.ac.il/PMN/
# format:
# tab-delimited with 4 columns <prot1> <prot2> <p-value> <direction>
# direction = 0 (undirected), 1 (directed forward), or 2 (directed backwards) 
# In pp interaction files, it will typically be 0 (undirected). 

##Set STRING confidence score filter
conf_filter = 900


#OHSU
setwd("~/Dropbox")

source("./thesis_work/code/thesis_funcs.R")

#Goto PPI data dir
setwd("./thesis_work/data/PPI/")

##########################################################
# 1. Import protein-protein interactions from two datasets
##########################################################

#Read in the STRING data, which are all functional associations for Mtb H37Rv from
#the STRING flatfile available on their website.

#Format is 3 columns, first two are the proteins and the last is the confidence score.
#Max confidence score is 1000.
string_raw <- read.table("./STRING/STRING.H37Rv.txt", stringsAsFactors=F)

#Original id looks like: 83332.Rv0002 so pull out just the Rv ID
string_raw$string1 <- sapply(string_raw[,1], function(x) strsplit(x, ".", fixed=T)[[1]][2])
string_raw$string2 <- sapply(string_raw[,2], function(x) strsplit(x, ".", fixed=T)[[1]][2])

string <- data.frame(string1=string_raw$string1, string2=string_raw$string2, conf=string_raw[,3],
                     stringsAsFactors=F)

#Read in the bacterial PPIs. First two rows are title and column labels, skip those.
#Only need first two columns which are the interacting proteins

hybrid <- read.table("pr100808n_si_008.txt", sep="\t", quote="", stringsAsFactors=F,
                     skip=2, colClasses=c("character", "character", rep("NULL",4)))

colnames(hybrid) <- c("hybrid1", "hybrid2")


##########################################################
# 2. Filter STRING associations for confidence > 0.9
##########################################################

#Filter for confidence score defined at top of script
string.filtered <- string[string$conf >= conf_filter, c(1,2)]

##Check how many genes overlap between data sets

hybrid.genes <- unique(c(hybrid$hybrid1, hybrid$hybrid2))
length(hybrid.genes)
#[1] 2907
string.filtered.genes <- unique(c(string.filtered$string1, string.filtered$string2))
length(string.filtered.genes)
#[1] 1678
length(intersect(hybrid.genes, string.filtered.genes))
#[1] 1203

string.700 <- string[string$conf >= 700, c(1,2)]
nrow(string.700)
#[1] 37094
string.700.genes <- unique(c(string.700$string1, string.700$string2))
length(string.700.genes)
#[1] 3301
length(intersect(hybrid.genes, string.700.genes))
#[1] 2424

##########################################################
# 3. Calculate the confidence score for each edge
##########################################################

#Need to create unique edge ID so we can see what edges
#are in common

string.filtered$edge_id <- paste(string.filtered$string1, string.filtered$string2, sep=",")
string.700$edge_id <- paste(string.700$string1, string.700$string2, sep=",")
hybrid$edge_id <- paste(hybrid$hybrid1, hybrid$hybrid2, sep=",")

#Check how many edges are found in both data sets
sum(string.filtered$edge_id %in% hybrid$edge_id)
#20 found in both

sum(string.700$edge_id %in% hybrid$edge_id)
#51

string.filtered[,1] <- as.character(string.filtered[,1])
string.filtered[,2] <- as.character(string.filtered[,2])

#At this point looks like adding a penalty for only being in one dataset
#will affect only a small part of the network, so not going to do that.

##USE 700 CONFIDENCE FILTER
colnames(hybrid)[1:2] <- c("e1", "e2")
colnames(string.700)[1:2] <- c("e1", "e2")

ppi.edges <- rbind(hybrid, string.700)
nrow(ppi.edges)
#[1] 45136
ppi.edges <- ppi.edges[!duplicated(ppi.edges$edge_id),c(1,2)]
nrow(ppi.edges)
#[1] 45085, took out 51 dupes

#Make PMN interaction file with weight based on confidence score

#Assign a confidence score of 900 to all 2-hybrid interactions
hybrid_conf <- hybrid[,c(1,2)]
hybrid_conf$conf <- 900

#Filter string dataset at 700 and keep confidence score
string.700_conf <- string[string$conf >= 700, ]
colnames(string.700_conf) <- c("e1", "e2", "conf")

ppi.edge_conf <- rbind(hybrid_conf, string.700_conf)
ppi.edge_conf$edge_id <- paste(ppi.edge_conf$e1, ppi.edge_conf$e2, sep=",")

nrow(ppi.edge_conf)
#[1] 45136
ppi.edge_conf <- ppi.edge_conf[!duplicated(ppi.edge_conf$edge_id), c(1,2,3)]

nrow(ppi.edge_conf)
# [1] 45085, took out 51 dupes

#Max confidence score of 1000 has weight of 0.05, 
#min confidence score of 700 has weight of 0.35
STRING_conf <- function(x){(1200 - x)/1000 - 0.15}

#Use combined degree from the two networks
node_ids <- unique(c(ppi.edges$e1, ppi.edges$e2))

length(node_ids)
#3784 genes

#initialize list to hold degrees
node_degrees <- vector("list", length(node_ids))
names(node_degrees) <- node_ids
node_degrees <- lapply(node_degrees, function(x) x <- 0)

#Add 1 to the degree of each protein in an edge
for (node in c(ppi.edges$e1, ppi.edges$e2)){
  node_degrees[[node]] <- as.integer(node_degrees[[node]]) + 1
}


#From paper, the confidence score for an edge is given as:
#
#  p = 1 / 2 + exp(-0.2x + 5)
#
#

calc_conf <- function(node_degree){
  p <- 1 / (2 + exp(-0.2*node_degree + 5))
  return(p)
}

calc_confs <- function(df, node_degrees){
  for (i in c(1:nrow(df))){
    e1.d <- node_degrees[[df[i,1]]]
    e2.d <- node_degrees[[df[i,2]]]
    degree <- max(e1.d, e2.d)
    df[i,3] <- calc_conf(degree)
  }
  
  return(df)
}

ppi.edges$conf <- 0

ppi.edges <- calc_confs(ppi.edges, node_degrees)

#Assign direction = 0, undirected (see PMN docs)
ppi.edges$direction <- 0

ppi.edges$e1 <- toupper(ppi.edges$e1)
ppi.edges$e2 <- toupper(ppi.edges$e2)

write.table(ppi.edges, file="mtb.pp.list", quote=F, sep="\t",
            col.names=F, row.names=F)
