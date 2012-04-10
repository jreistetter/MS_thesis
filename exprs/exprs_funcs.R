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
    probe_avg <- unlist(lapply(probes[,1:(ncol(df)-1)], function(x) mean(x, na.rm=T)))
    avg_df[i,1:(ncol(df)-1)] <- probe_avg
    avg_df[i,]$gene <- gene
    
    i = i+1
  }
  colnames(avg_df) <- colnames(df)
  final_df <- avg_df[,1:(ncol(df)-1)]
  rownames(final_df) <- avg_df$gene
  return(final_df)
  
}

discretizer <- function(val, fold){
  if(is.na(val)){
    return(NA)
  }
  if(val >= fold){
    return(1)
  }
  
  if(val <= -1*fold){
    return (-1)
  }
  return(0)
}

vec.discret <- function(vec, fold){
  return(vapply(vec, discretizer, FUN.VALUE=c(1),
                USE.NAMES=F, fold = fold))
}

discretize <- function(df, fold=1.5){
  disc.df <- data.frame(lapply(df, vec.discret, fold=fold))
  rownames(disc.df) <- rownames(df)
  return(disc.df)
}

probe_mean_A <- function(gene_id, df){
  #get intensities for all arrays of that gene
  probes.A <- df[df$gene == gene_id,1:(ncol(df)-1)]
  mean.A <- apply(probes.A, 2, mean)
  return(mean.A)
}


probe_CV <- function(gene_id, df){
  #get intensities for all arrays of that gene
  probes.M <- df[df$gene == gene_id,1:(ncol(df)-1)]
  coef_vars <- apply(probes.M, 2, coef_var)
  return(coef_vars)
}

coef_var <- function(vals){
  coef <- sd(vals, na.rm=T) / mean(vals, na.rm=T)
  return(coef)
}

remove_bad_spots <- function(ma_list, type="M"){
  probe.weights <- ma_list$weights
  probe.weights[probe.weights == 0] <- NA
  if (type=="M"){
    cleaned <- probe.weights * ma_list$M
  }
  
  if (type=="A"){
    cleaned <- probe.weights * ma_list$A
  }
  rownames(cleaned) <- ma_list$genes$Name
  return(as.data.frame(cleaned))

}

get_consensus <- function(vals, threshold){
  #Takes a vector of log2(R/G) values and finds the
  #consensus, defined as nvals - 1 having the same
  #call.
  #
  #Params:
  #vals - vector of log2(R/G) values
  #threshold - threshold for calling up/down/nocall for expression
  n_vals <- sum(!is.na(vals))
  if (sum(vals >= threshold) >= n_vals-1){
    return(1)
  }
  
  if(sum(vals <= -1*threshold) >= n_vals-1){
    return(-1)
  }
  if(sum(vals >= -1*threshold & vals <= threshold) >= n_vals-1){
    return(0)
  }
  return(NA)
}

#Test the function
stopifnot(get_consensus(c(0.1, 0.2, 0, -0.1), 1) == 0)
stopifnot(get_consensus(c(1.1, 1.2, 0.9, 1.3), 1) == 1)
stopifnot(get_consensus(c(-1.1, -1.2, -1.8, -0.1), 1) == -1)
stopifnot(is.na(get_consensus(c(-1.1, 0.2, 1.3, 1.8), 1)))
stopifnot(is.na(get_consensus(c(-0.9, 0.2, 1.3, 1.8), 1)))





