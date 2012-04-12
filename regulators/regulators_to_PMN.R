# Script to parse files downloaded from AmiGO for various GO categories
#
# Written by Joe Reistetter

# GO categories:
#   GO:0006950 : response to stress [188 gene products]
#   GO:0006955 : immune response [25 gene products]
#   GO:0009607 : response to biotic stimulus [122 gene products]
#   GO:0042221 : response to chemical stimulus [140 gene products]
#   GO:0051716 : cellular response to stimulus [96 gene products]
#
#   GO:0005576 : extracellular region [276 gene products] 
#   GO:0044425 : membrane part [104 gene products]

#laptop
setwd("~/schoolDB/Dropbox/thesis_work/data/regulators/")


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
  return(g.rv)
}
 

#parse in data downloaded from AmiGO
g.6950.rv <- parse_GO("GO_0006950.txt")
g.6955.rv <-  parse_GO("GO_0006955.txt")
g.9607.rv <- parse_GO("GO_0009607.txt")
g.42221.rv <- parse_GO("GO_0042221.txt")
g.51716.rv <- parse_GO("GO_0051716.txt")

g.5576.rv <- parse_GO("GO_0005576.txt")
g.44425.rv <- parse_GO("GO_0044425.txt")

go.rv <- unique(c(g.42221.rv, g.44425.rv, g.51716.rv,
                  g.5576.rv, g.6950.rv, g.6955.rv,
                  g.9607.rv))

go.rv <- go.rv[!is.na(go.rv)]

length(go.rv)
#[1] 622

length(c(g.42221.rv, g.44425.rv, g.51716.rv,
         g.5576.rv, g.6950.rv, g.6955.rv,
         g.9607.rv))
#[1] 951
save(go.rv, file="go.rv.RData")
write.table(go.rv, file="go_regulators.txt", quote=F,
row.names=F, col.names=F)

#Combine with other regulator group

reg_list <- read.table("initial.regulators", head=F, stringsAsFactors=F)

go.rv <- toupper(go.rv)

final_reg_list <- unique(c(go.rv, toupper(reg_list[,1])))

length(final_reg_list)
#[1] 650

