#Script to load the parsed PMN objects into R for analysis
#Data was parsed from PMN output in 4_17_30_mods_parse.py

#Written by Joe Reistetter

options(stringsAsFactors=F)
#load data
setwd("~/Dropbox/thesis_work/PMN_output/")

cpds <- read.table("4.17.30_mods_parsed.txt",
                   head=T, sep='\t', stringsAsFactors=F)

#Extract list of parents
comb_par <- unique(cpds$parents)
parents <- unlist(sapply(comb_par, function(x) strsplit(x, " ", fixed=T), USE.NAMES=F))
parents <- parents[!is.na(parents)]

#Module membership
mod_members <- read.table("4.17.30_mods_members.txt",
                       head=T, sep='\t', stringsAsFactors=F)

#Number of genes per module
mod_sizes <- as.data.frame(table(mod_members$moduleID))
colnames(mod_sizes) <- c("moduleID", "n_genes")

#Pathways
path_df <- read.table("4.17.30_mods_pathways.txt",
                      head=F, sep=':', stringsAsFactors=F)

#Parse into a list
pathways <- list()

for (i in 1:nrow(path_df)){
  mod_id = path_df[i,1]
  pathways[[mod_id]] <- unlist(strsplit(path_df[i,2], "\t", fixed=T))
}

#Calculate the size of each pathway
sizes <- unlist(lapply(path_sizes, length))
ids <- names(sizes)

pathway_sizes <- data.frame(moduleID=ids, path_size=sizes)

#Merge pathway size and module size together
stats <- merge(mod_sizes, pathway_sizes)


#Write out data
write.table(stats, file="4.17_30mods_genes_pathsizes.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

write.table(cpds, file="4.17_30mods_CPDs.txt",
            col.names=T, row.names=F, quote=F, sep='\t')
