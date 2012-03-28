# Script containing support functions for the parsing, processing, and writing
# of the pDNA interaction data.


assign_pDNA_edge <- function(reg_target, genes, op_list){
  #Assigns edges between regulator and targets, and members of operons' targets
  #that are downstream of target
  #Params:
#   reg_target - df, 1st col is regulator, 2nd target
#   genes - df, 1st col is operon ID corresponding to list name, 2nd is gene
#   op_list - list of vectors with operon members, key is the operon ID from genes
  
  pDNA_edges <- as.matrix(data.frame(protein=vector(mode="character"),
                                     gene=vector(mode="character"),
                                     stringsAsFactors=F))
  
  for (i in 1:nrow(reg_target)){
    target <- reg_target[i,2]
    
    #If target isn't in an operon, just write that pair
    if(!(target%in%genes[,2])){
      pDNA_edges <- rbind(pDNA_edges, c(reg_target[i,1], target))
      next
    }
    
    op.id <- genes[genes[,2]==target,1]
    op.vec <- op_list[[op.id]]
    #Get index of target in operon
    target.idx <- which(toupper(op.vec)==target)
    #Check if on complement strand
    complement = grepl("C", target, ignore.case=TRUE)
    
    #If on c strand, need to add target and all 
    #genes to the "left" of the target
    if (complement){
      for (j in 1:target.idx){
        pDNA_edges <- rbind(pDNA_edges, c(reg_target[i,1], op.vec[j]))
      }
      next
    }
    
    #If not on c strand, add target and all genes to "right" of target
    for (k in target.idx:length(op.vec)){
      pDNA_edges <- rbind(pDNA_edges, c(reg_target[i,1], op.vec[k]))
    }
  }
  return(pDNA_edges)
}

operons_from_gene_pairs <- function(op.df, gene.cols){
  #Params
  #op.df - data frame with operon gene pairs as two columns
  #gene.cols - character vector of the column indices of the gene IDs
  
  #Initializes variables to use in loop that creates operons
  op.list <- list()
  i <- 1
  op_id <- 1
  
  #Loop through each row of the filtered operon predictions.
  #Depending on the continuity of the operon, loop does
  #different things.
  while (i < nrow(op.df)+1){
    
    #If there is an NA, skip that row
    if(is.na(op.df[i+1,gene.cols[1]]) | is.na(op.df[i,gene.cols[2]])){
      i <- i+1
      op_id <- op_id+1
      next
    }
    
    #If the operon is only 2 genes, add those 2 genes as a new operon
    if(op.df[i+1,gene.cols[1]] != op.df[i,gene.cols[2]]){
      op_id <- op_id+1
      op.list[[as.character(op_id)]] <- c(op.list[[as.character(op_id)]],
                                               unlist(c(op.df[i,gene.cols])))
      i <- i+1
      op_id <- op_id+1
      next
    }
    
    #While the genes are continuous (in an operon), keep addding to new operon
    while(op.df[i+1,gene.cols[1]] == op.df[i,gene.cols[2]]){
      
      existing <- op.list[[as.character(op_id)]]
      these_genes <- unlist(c(op.df[i,gene.cols], op.df[i+1,gene.cols]))
      genes_add <- c(existing, these_genes)
      op.list[[as.character(op_id)]] <- genes_add
      i <- i+1
    }
    op_id <- op_id+1
    i <- i+1
  }
  
  #Add last operon
  op.list[[length(op.list)+1]] <- 
    unlist(c(op.df[i-1,gene.cols]))
  
  #For simplicity, the while loop doesn't check if a gene is already in an operon
  #before adding it. This removes duplicates. Based on algorithm, duplicates
  #are expected..,.,.. 
  
  op.list <- lapply(op.list, unique)
  
  return(op.list)
}

op_list_to_df <- function(op.list){
  #Takes a list of operons as produced by
  #operons_from_gene_pairs and returns a dataframe
  #with each gene and its operon ID.
  #Params
  #op.list - list of operons produced by operons_from_gene_pairs

  df <- data.frame(operon_ID=c(), gene=c())
  
  for (i in c(1:length(op.list))){
    genes <- op.list[[i]]
    rows <- matrix(c(rep(i, length(genes)), genes), ncol=2)
    df <- rbind(df, rows)
  }
  
  colnames(df) <- c("operon_ID", "gene_ID")
  
  return(df)
  }