#Functions used across scripts to process expression data

# avg_probes -- Average intensities for arrays with multiple probes per gene
# vec_discret -- vectorized version of discretize to use on a dataframe of intensities
# discretize -- discretizes a vector of intensities to {-1,0,1}

avg_probes <- function(df, gene_ids){
  #initialize datafame to hold results
  avg_df <- as.data.frame(
    matrix(nrow=length(gene_ids), ncol=ncol(df)-1) 
    )
  avg_df$gene <- ""
  
  i = 1
  for (gene in gene_ids){
    probes <- df[which(df$gene == gene),]
    probe_avg <- unlist(lapply(probes[,1:(ncol(df)-1)], function(x){sum(x)/length(x)}))
    avg_df[i,1:(ncol(df)-1)] <- probe_avg
    avg_df[i,]$gene <- gene
    
    i = i+1
  }
  colnames(avg_df) <- colnames(df)
  final_df <- avg_df[,1:(ncol(df)-1)]
  rownames(final_df) <- avg_df$gene
  return(final_df)
  
}

discretizer <- function(val){
  if(val >= 1){
    return(1)
  }
  
  if(val <= -1){
    return (-1)
  }
  return(0)
}

vec.discret <- function(vec){
  return(sapply(vec, discretizer, USE.NAMES=F))
}

discretize <- function(df){
  data.frame(lapply(df, vec.discret))
}