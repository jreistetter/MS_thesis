#Script to load the parsed PMN objects into R for analysis
#Data was parsed from PMN output in 8.10.12_30_mods_0.1_highvar_parse.py

#Written by Joe Reistetter

options(stringsAsFactors=F)
#load data
setwd("~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/")

cpds <- read.table("8.10.12_30_mods_0.1_highvar_mods_parsed.txt",
                   head=T, sep='\t', stringsAsFactors=F)

#Find modules that have all uniform CPDs
non_uniform <- c()
for (mod in unique(cpds$moduleID)){
  mod_cpds <- unlist(cpds[cpds$moduleID == mod, 4:6])
  if (!all(mod_cpds == 0.333333, na.rm=T)){
    non_uniform <- c(non_uniform, mod)
  }
}

length(non_uniform)
#[1] 28, only removed 2 modules that way.
all_mod_ids <- unique(cpds$moduleID)
all_mod_ids[which(!(all_mod_ids %in% non_uniform))]
#[1] "mod29 " "mod10 "

cpds.all_not_uniform <- cpds[cpds$moduleID %in% non_uniform,]
cpds.all_not_uniform$uniform <- TRUE
#Filter out to get only the non-uniform probs
for (i in 1:nrow(cpds.all_not_uniform)){
  cpds.all_not_uniform[i,]$uniform <- all(cpds.all_not_uniform[i,4:6]==0.333333)
}

cpds.non_uniform <- cpds.all_not_uniform[cpds.all_not_uniform$uniform==F,]
mod_n_non_uni <- as.data.frame(table(cpds.non_uniform$moduleID))
colnames(mod_n_non_uni) <- c("moduleID", "n_non_uni")

count_changed <- function(cpds, threshold){
  mod_n_up_or_down <- list()
  
  for (mod in unique(cpds$moduleID)){
    mod_cpds <- cpds[cpds$moduleID==mod,]
    n_down <- sum(mod_cpds$mod_down > threshold)
    n_up <- sum(mod_cpds$mod_up > threshold)
    mod_n_up_or_down[[mod]] <- sum(n_down, n_up)
  }
  
  mod_n_dif <- data.frame(moduleID=names(mod_n_up_or_down),
                          n_changed=unlist(mod_n_up_or_down),
                          stringsAsFactors=F)
  colnames(mod_n_dif)[2] <- paste("thresh.", threshold, sep="")
  
  return(mod_n_dif)
}

n_changed.2 <- count_changed(cpds.non_uniform, 0.2)
n_changed.3 <- count_changed(cpds.non_uniform, 0.3)
n_changed.4 <- count_changed(cpds.non_uniform, 0.4)
n_changed.5 <- count_changed(cpds.non_uniform, 0.5)
n_changed.6 <- count_changed(cpds.non_uniform, 0.6)
n_changed.7 <- count_changed(cpds.non_uniform, 0.7)
n_changed.8 <- count_changed(cpds.non_uniform, 0.8)




#Module membership
mod_members <- read.table("8.10.12_30_mods_0.1_highvar_mods_members.txt",
                       head=T, sep='\t', stringsAsFactors=F)

#Number of genes per module
mod_sizes <- as.data.frame(table(mod_members$moduleID))
colnames(mod_sizes) <- c("moduleID", "n_genes")

#Pathways
path_df <- read.table("8.10.12_30_mods_0.1_highvar_mods_pathways.txt",
                      head=F, sep=':', stringsAsFactors=F)

#Parse into a list
pathways <- list()

for (i in 1:nrow(path_df)){
  mod_id = path_df[i,1]
  pathways[[mod_id]] <- unlist(strsplit(path_df[i,2], "\t", fixed=T))
}

#Calculate the size of each pathway
sizes <- unlist(lapply(pathways, length))
ids <- names(sizes)

pathway_sizes <- data.frame(moduleID=ids, path_size=sizes)

#Merge pathway size and module size, and n_changed together
stats <- merge(mod_sizes, pathway_sizes)
stats <- merge(stats, mod_n_non_uni)
stats <- merge(stats, n_changed.2)
stats <- merge(stats, n_changed.3)
stats <- merge(stats, n_changed.4)
stats <- merge(stats, n_changed.5)
stats <- merge(stats, n_changed.6)
stats <- merge(stats, n_changed.7)
stats <- merge(stats, n_changed.8)

#Write out data
write.table(stats, 
            file="./output/8.10.12_30_mods_0.1_highvar_genes_pathsizes.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

write.table(cpds, file="./output/8.10.12_30_mods_0.1_highvar_CPDs.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

#Write out a file for each module for submission to GOEAST, including parents
comb_par <- unique(cpds.non_uniform$parents)
parents <- unlist(sapply(comb_par, function(x) strsplit(x, " ", fixed=T), USE.NAMES=F))
parents <- parents[!is.na(parents)]

for (mod_id in non_uniform){
  mod_genes <- mod_members[mod_members$moduleID==mod_id,2]
  mod_parents_raw <- cpds.non_uniform[cpds.non_uniform$moduleID==mod_id,2][1]
  mod_parents <- unlist(strsplit(mod_parents_raw, " ", fixed=T))
  mod_genes <- c(mod_parents, mod_genes)
  f_name <- paste("./output/GOEAST/", mod_id, "_genes.txt", sep="")
  write.table(mod_genes, f_name, quote=F, row.names=F, col.names=F)
  
}

