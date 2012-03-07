# Script containing support functions for the parsing, processing, and writing
# of the pDNA interaction data.

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