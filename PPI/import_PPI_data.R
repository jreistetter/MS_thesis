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
#setwd("~/Dropbox")

#Laptop
setwd("~/schoolDB/Dropbox")

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

string <- data.frame(string1=string_raw$string1, string2=string_raw$string2, conf=string_raw[,3])

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

##########################################################
# 3. Calculate the confidence score for each edge
##########################################################

#Need to create unique edge ID so we can see what edges
#are in common

string.filtered$edge_id <- paste(string.filtered$string1, string.filtered$string2, sep=",")
hybrid$edge_id <- paste(hybrid$hybrid1, hybrid$hybrid2, sep=",")

#Check how many edges are found in both data sets
sum(string.filtered$edge_id %in% hybrid$edge_id)
#20 found in both

string.filtered[,1] <- as.character(string.filtered[,1])
string.filtered[,2] <- as.character(string.filtered[,2])

#At this point looks like adding a penalty for only being in one dataset
#will affect only a small part of the network, so not going to do that.

#Use combined degree from the two networks
node_ids <- unique(c(string.filtered$string1, string.filtered$string2, 
                     hybrid$hybrid1, hybrid$hybrid2))

length(node_ids)
#3382 genes

#initialize list to hold degrees
node_degrees <- vector("list", length(node_ids))
names(node_degrees) <- node_ids
node_degrees <- lapply(node_degrees, function(x) x <- 0)

#Add 1 to the degree of each protein in an edge for STRING
for (node in c(string.filtered$string1, string.filtered$string2)){
  node_degrees[[node]] <- as.integer(node_degrees[[node]]) + 1
}

#Add 1 to the degree of each protein in an 2-hybrid edge
for (node in c(hybrid$hybrid1, hybrid$hybrid2)){
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

string.filtered$conf <- 0
add_degrees <- function(df, node_degrees){
  for (i in c(1:nrow(df))){
    e1.d <- node_degrees[[df[i,1]]]
    e2.d <- node_degrees[[df[i,2]]]
    degree <- max(e1.d, e2.d)
    df[i,4] <- calc_conf(degree)
  }
  
  return(df)
}

string.conf <- add_degrees(string.filtered, node_degrees)
hybrid.conf <- add_degrees(hybrid, node_degrees)

colnames(hybrid.conf) <- c("e1", "e2", "edge_id", "conf")
colnames(string.conf) <- c("e1", "e2", "edge_id", "conf")

ppi.edges.dupes <- rbind(hybrid.conf[,c(1,2,4)], string.conf[,c(1,2,4)])

dupes <- duplicated(ppi.edges.dupes)

ppi.edges <- ppi.edges.dupes[!dupes,]

nrow(ppi.edges.dupes) - nrow(ppi.edges)
#Should equal 20, the number of edges in common between
#the two networks

#Assign direction = 0, undirected (see PMN docs)
ppi.edges$direction <- 0
write.table(ppi.edges, file="mtb.pp.list", quote=F, sep="\t",
            col.names=F, row.names=F)
