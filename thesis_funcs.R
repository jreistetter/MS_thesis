gene.mod <- function(geneID, modules){
  return(modules[modules$gene==geneID,]$moduleID)
}

get_module <- function(modID, modules){
  #Returns char vec of genes belonging to modID
  return(modules[modules$moduleID==modID,2])
}

get_tf <- function(modID, modules, pDNA){
  mod.members <- get_module(modID, modules)
  mod.tfs <- table(pDNA[pDNA$target%in%mod.members,]$tf)
  mod.tfs.df <- as.data.frame(mod.tfs)
  if (nrow(mod.tfs.df)==0){
    return(data.frame(Var1=c(NA), Freq=c(NA), modID=c(modID)))
  }
  mod.tfs.df$moduleID <- modID
  return(mod.tfs.df)
}