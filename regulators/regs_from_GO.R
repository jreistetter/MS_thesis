# Script to parse files downloaded from AmiGO for various GO categories
#
# Written by Joe Reistetter

# GO categories:
#   GO:0006950 : response to stress [188 gene products]

#laptop
setwd("~/schoolDB/Dropbox/thesis_work/data/regulators/")


get_rv <- function(pipe_list){
  #Extracts Rv ID from pipe-delimted list
  ids <- unlist(strsplit(pipe_list, '|', fixed=T))
  rv_id <- ids[grepl("Rv", ids, fixed=T)]
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
  return(g.rv)
}

###################################
#
#     GO:0006950 : response to stress [188 gene products]
#
###################################

g.6950.rv <- parse_GO("GO_0006950.txt")
