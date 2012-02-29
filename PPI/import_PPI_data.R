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

#Read in the STRING data, which are all functional associations for Mtb H37Rv from
#the STRING flatfile available on their website.

#Format is 3 columns, first two are the proteins and the last is the confidence score.
#Max confidence score is 1000.


888.Rv

string_raw <- read.table("./STRING/STRING.H37Rv.txt", stringsAsFactors=F)

#Original id looks like: 83332.Rv0002 so pull out just the Rv ID
string_raw$string1 <- sapply(string_raw[,1], function(x) strsplit(x, ".", fixed=T)[[1]][2])
string_raw$string2 <- sapply(string_raw[,2], function(x) strsplit(x, ".", fixed=T)[[1]][2])

string <- data.frame(string1=string_raw$string1, string2=string_raw$string2, conf=string_raw[,3])

#Filter for confidence score defined at top of script
string.filtered <- string[string$conf >= conf_filter, c(1,2)]

#Read in the bacterial PPIs. First two rows are title and column labels, skip those.
#Only need first two columns which are the interacting proteins

hybrid <- read.table("pr100808n_si_008.txt", sep="\t", quote="", stringsAsfactors=F
                     skip=2, colClasses=c("character", "character", rep("NULL",4)))

colnames(hybrid) <- c("hybrid1", "hybrid2")


