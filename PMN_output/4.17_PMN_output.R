#Script to load the parsed PMN objects into R for analysis
#Data was parsed from PMN output in 4_17_30_mods_parse.py

#Written by Joe Reistetter

options(stringsAsFactors=F)
#load data
setwd("~/Dropbox/thesis_work/PMN_output/")

cpds <- read.table("4.17.30_mods_parsed.txt",
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
#[1] "mod30 " "mod26 "

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
mod_members <- read.table("4.17.30_mods_members.txt",
                       head=T, sep='\t', strngsAsFactors=F)

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
write.table(stats, file="4.17_30mods_genes_pathsizes.txt",
            col.names=T, row.names=F, quote=F, sep='\t')

write.table(cpds, file="4.17_30mods_CPDs.txt",
            col.names=T, row.names=F, quote=F, sep='\t')


#Extract list of parents
comb_par <- unique(cpds$parents)
parents <- unlist(sapply(comb_par, function(x) strsplit(x, " ", fixed=T), USE.NAMES=F))
parents <- parents[!is.na(parents)]

#Extract list of parents of modules with 1 up/down CPD
good_mods <- c("mod5", "mod28", "mod23", "mod25", "mod18", "mod19")
good_par <- unlist(
  sapply(unique(cpds[cpds$moduleID %in% good_mods,]$parents), 
         function(x) strsplit(x, " ", fixed=T),
         USE.NAMES=F))

#See if those parents are in other verified modules
sum(good_par %in% mod_members[mod_members$moduleID %in% good_mods,]$gene)
#[1] 0, so no regulators of a module are in another verified module


#Read in regulator lists to see where they came from
get_rv <- function(pipe_list){
  #Extracts Rv ID from pipe-delimted list
  ids <- unlist(strsplit(pipe_list, '|', fixed=T))
  rv_id <- ids[grepl("Rv", ids, fixed=T)]
  if (length(rv_id) == 0){
    return(NA)
  }
  if (length(rv_id) > 1){
    return(rv_id[2])
  }
  return(rv_id)
}

parse_GO <- function(go_path){
  g.raw <- read.table(go_path, head=F, stringsAsFactors=F,
                      quote="", sep="\t")
  
  #Pull out two rows with name and pipe-delimited list of different gene IDs
  g.cols <- g.raw[,c(3,11)]
  colnames(g.cols) <- c("name", "id_list")
  
  g.rv <- unique(sapply(g.cols$id_list, get_rv))
  return(toupper(g.rv))
}

setwd("~/Dropbox/thesis_work/data/regulators/")
init.reg <- as.character(read.table("initial.regulators", head=F,
                       stringsAsFactors=F)[,1])

all.reg <- as.character(read.table("final.regulators",
                      head=F, stringsAsFactors=F)[,1])

g.6950.rv <- parse_GO("GO_0006950.txt")
g.6955.rv <-  parse_GO("GO_0006955.txt")
g.9607.rv <- parse_GO("GO_0009607.txt")
g.42221.rv <- parse_GO("GO_0042221.txt")
g.51716.rv <- parse_GO("GO_0051716.txt")

g.5576.rv <- parse_GO("GO_0005576.txt")
g.44425.rv <- parse_GO("GO_0044425.txt")

parents_hyper <- function(parents, regs, all_regs){
  n_parents <- length(parents)
  n_not <- length(all_regs) - n_parents
  reg_parents <- sum(parents %in% regs)
  n_regs_not <- length(regs) - reg_parents
  pval <- phyper(reg_parents, n_parents, n_not, n_regs_not)
  return(1-pval)
}

n_par <- length(parents)
n_not <- length(all.reg) - length(parents)

sum(parents%in%init.reg)
sum(parents%in%g.6950.rv)
sum(parents%in%g.6955.rv)
sum(parents%in%g.9607.rv)
sum(parents%in%g.42221.rv) 
sum(parents%in%g.51716.rv)
sum(parents%in%g.5576.rv)
sum(parents%in%g.44425.rv)

parents_hyper(parents, init.reg, all.reg)
parents_hyper(parents, g.6950.rv, all.reg)
parents_hyper(parents, g.9607.rv, all.reg)
parents_hyper(parents, g.42221.rv, all.reg)
parents_hyper(parents, g.51716.rv, all.reg)
parents_hyper(parents, g.5576.rv, all.reg)
parents_hyper(parents, g.44425.rv, all.reg)





