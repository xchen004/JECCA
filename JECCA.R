JECCA <- function(da,metadata,ncomp=NULL)
{
  # This function perform CCA analysis for multiple datasets
  #
  # Args:
  #   da:       gene expression matrix of combined multiple datasets.
  #             Each column represents a sample, each row represents a gene..
  #   metadata: meta file with first column of sample id matched with 
  #             the samples in gene expression matrix and the second 
  #             column indicating the dataset information for each sample.
  #   ncomp:    total number of JECCA embeddings need to be returned
  #             The number of embeddings should not exceed the minimum number
  #             of samples of the multiple datasets.
  # Returns:
  #   res:      JECCA derived embeddings
  #           
  
  require(PMA)
  rownames(metadata) <- metadata$sid
  sid <- intersect(colnames(da), metadata$sid)
  if(length(sid) == 0)
  {stop('Not matching samples between data matrix and meta file')}
  if(length(sid) > 0)
  {
    ds <- unique(metadata$DataSet)
    da <- da[,sid]
    metadata <- metadata[sid,]
    da_s <- sn <- NULL
    for(i in 1:length(ds))
    {
      temp_anno <- metadata[metadata$DataSet %in% ds[i],]
      temp_da <- as.matrix(da[,rownames(temp_anno)])
      da_s[[i]] <- t(scale(t(temp_da)))
      da_s[[i]][is.nan(da_s[[i]])] <- 0
      sn <- c(sn, ncol(temp_da))
    }
    if(is.null(ncomp))
    {ncomp <- min(sn)}
    if(ncomp > min(sn)){
      ncomp <- min(sn)
      print('Exceeding permissible dimensions, employing the maximal ncomp as a substitute')
    }
    set.seed(123)
    perm.out <- PMA::MultiCCA.permute(da_s)
    out <- PMA::MultiCCA(da_s, penalty=perm.out$penalties[,10], 
                         ncomponents=ncomp, ws=perm.out$ws.init)
    for(i in 1:length(da_s))
    {
      rownames(out$ws[[i]]) <- colnames(da_s[[i]])
      colnames(out$ws[[i]]) <- paste0('JECCA',1:ncomp)
    }
    res <- do.call(rbind, out$ws)
    return(res)
  }
}




