JECCA_pathway <- function(da,metadata,embedding,pathway_database,output)
{
  # This function perform CCA analysis for multiple datasets
  #
  # Args:
  #   da:               gene expression matrix of combined multiple datasets.
  #                     Each column represents a sample, each row represents a gene.
  #   metadata:         meta file with first column of sample id matched with 
  #                     column names of gene expression matrix and the second column 
  #                     indicating the dataset information for each sample.
  #   embedding:        JECCA output embeddings based on da matrix.            
  #   pathway_database: The pathway database need to be selected for the GSEA pathways analysis. 
  #                     The following pathway databases are provided: ReactomPA, KEGG, GOALL, GOBP, GOMF, GOALL. 
  #                     User can select any combinations of the above databases to perform GSEA pathway analysis. 
  #                     fgsea function is used to perform reactomPA GSEA analysis. 
  #                     gseKEGG and gseGO function is used to perform KEGG and GO GSEA analysis.
  #   output:           The output directory. The function will generate following
  #                     results in the output directory: subtype expression based on the da matrix and 
  #                     embedding, GSEA pathway results.
  #           
  require(biomaRt)
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2)
  require(AWFisher)
  require(dorothea)
  require(reshape2)
  require(fgsea)
  require(org.Hs.eg.db)
  entrez_to_symbol <- function(gene, anno)
  {
    gene <- gene
    anno <- as.data.frame(anno)
    entrz <- unlist(strsplit(gene, '/'))
    symbol <- anno[match(entrz, anno$entrez),]$gene
    return(paste(symbol,collapse='/'))
  }
  
  rownames(metadata) <- metadata$sid
  sid <- intersect(colnames(da), metadata$sid)
  if(length(sid) == 0)
  {stop('Not matching samples between data matrix and meta file')}
  if(length(sid) > 0)
  {
    ds <- unique(metadata$DataSet)
    da <- da[,sid]
    metadata <- metadata[sid,]
    da_s <- NULL
    for(i in 1:length(ds))
    {
      temp_anno <- metadata[metadata$DataSet %in% ds[i],]
      temp_da <- as.matrix(da[,rownames(temp_anno)])
      temp_s <- t(scale(t(temp_da)))
      temp_s[is.nan(temp_s[[i]])] <- 0
      da_s <- cbind(da_s,temp_s)
    }
    da_s <- scale(da_s)
    ## generate subtype expression profile 
    subtype_exprs <- c()
    jecca_names <- paste0('JECCA',1:ncol(embedding))
    for(i in 1:length(ds))
    {
      temp_exprs <- as.matrix(da_s[,rownames(metadata[metadata$DataSet == ds[i],])]) %*%
        as.matrix(embedding[rownames(metadata[metadata$DataSet == ds[i],]),])
      colnames(temp_exprs) <- paste0(ds[i],'_',jecca_names)
      if(i == 1)   subtype_exprs <- temp_exprs
      if(i > 1) subtype_exprs <- cbind.data.frame(subtype_exprs, temp_exprs)
    }
    subtype_exprs2 <- c()
    for(i in 1:length(jecca_names))
    {
      index <- seq(i,ncol(subtype_exprs), by=length(jecca_names))
      if(i == 1) subtype_exprs2 <- subtype_exprs[,index]
      if(i > 1) subtype_exprs2 <- cbind.data.frame(subtype_exprs2,subtype_exprs[,index])
    }
    subtype_exprs <- subtype_exprs2
    hg_genes <- rownames(subtype_exprs)
    subtype_exprs <- cbind.data.frame(gene=hg_genes, subtype_exprs)
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    geneid <- getBM(attributes = c('entrezgene_id', 'external_gene_name'),
                    filters = 'external_gene_name', values = hg_genes, mart = ensembl)
    geneid <- geneid[!is.na(geneid[,1]),]
    geneid <- geneid[!duplicated(geneid[,1]),]
    colnames(geneid) <- c('entrez', 'gene')
    subtype_exprs <- merge(geneid, subtype_exprs, by = 'gene')
    subtype_exprs <- subtype_exprs[!duplicated(subtype_exprs$entrez),]
    write.csv(subtype_exprs,file=paste0(output,'subtype_exprs.csv'),row.names=F)
    
    ## pathway analysis ##
    for(i in pathway_database)
    {assign(paste0("gse_",i), NULL)}
    for(i in 3:ncol(subtype_exprs))
    {
      gene_list <- subtype_exprs[,i]
      # name the vector
      names(gene_list) <- subtype_exprs$entrez
      # omit any NA values 
      gene_list<-na.omit(gene_list)
      # sort the list in decreasing order (required for clusterProfiler)
      gene_list = sort(gene_list, decreasing = TRUE)
      
      if('ReactomPA' %in% pathway_database)
      {
        ## get reactomPA database
        ra_pathways <- reactomePathways(as.character(subtype_exprs$entrez))
        gse_ReactomPA[[i-2]] <- fgsea(ra_pathways, gene_list, minSize=10, maxSize=500)
        gene_symbol <- c()
        for(j in 1:nrow(gse_ReactomPA[[i-2]]))
        {
          gse_ReactomPA[[i-2]]$leadingEdge[j] <- paste(unlist(gse_ReactomPA[[i-2]]$leadingEdge[j]),
                                                       collapse='/')
          gene_symbol <- c(gene_symbol, entrez_to_symbol(unlist(gse_ReactomPA[[i-2]][j,8]), geneid))
        }
        gse_ReactomPA[[i-2]] <- cbind.data.frame(gse_ReactomPA[[i-2]], GeneID2=gene_symbol)
      }
 
      if('KEGG' %in% pathway_database)
      {
        gse_KEGG[[i-2]] <- gseKEGG(geneList=gene_list, organism='hsa',
                                   pvalueCutoff = 1, verbose = FALSE,
                                   pAdjustMethod = "fdr")
        gene_symbol <- c()
        for(j in 1:nrow(gse_KEGG[[i-2]]))
        {
          gene_symbol <- c(gene_symbol,entrez_to_symbol(gse_KEGG[[i-2]][j,11], geneid))
        }
        gse_KEGG[[i-2]] <- cbind.data.frame(gse_KEGG[[i-2]], GeneID2=gene_symbol)
      }

      if('GOALL' %in% pathway_database)
      {
        gse_GOALL[[i-2]] <- gseGO(geneList=gene_list, ont='ALL',
                                  OrgDb=org.Hs.eg.db, keyType='ENTREZID',
                                  pvalueCutoff = 1, verbose = FALSE,
                                  pAdjustMethod = "fdr")
        gene_symbol <- c()
        for(j in 1:nrow(gse_GOALL[[i-2]]))
        {
          gene_symbol <- c(gene_symbol,entrez_to_symbol(gse_GOALL[[i-2]][j,11], geneid))
        }
        gse_GOALL[[i-2]] <- cbind.data.frame(gse_GOALL[[i-2]], GeneID2=gene_symbol)
      }
 
      if('GOBP' %in% pathway_database)
      {
        gse_GOBP[[i-2]] <- gseGO(geneList=gene_list, ont='BP',
                                  OrgDb=org.Hs.eg.db, keyType='ENTREZID',
                                  pvalueCutoff = 1, verbose = FALSE,
                                  pAdjustMethod = "fdr")
        gene_symbol <- c()
        for(j in 1:nrow(gse_GOBP[[i-2]]))
        {
          gene_symbol <- c(gene_symbol,entrez_to_symbol(gse_GOBP[[i-2]][j,11], geneid))
        }
        gse_GOBP[[i-2]] <- cbind.data.frame(gse_GOBP[[i-2]], GeneID2=gene_symbol)
      }
      
      if('GOMF' %in% pathway_database)
      {
        gse_GOMF[[i-2]] <- gseGO(geneList=gene_list, ont='MF',
                                 OrgDb=org.Hs.eg.db, keyType='ENTREZID',
                                 pvalueCutoff = 1, verbose = FALSE,
                                 pAdjustMethod = "fdr")
        gene_symbol <- c()
        for(j in 1:nrow(gse_GOMF[[i-2]]))
        {
          gene_symbol <- c(gene_symbol,entrez_to_symbol(gse_GOMF[[i-2]][j,11], geneid))
        }
        gse_GOMF[[i-2]] <- cbind.data.frame(gse_GOMF[[i-2]], GeneID2=gene_symbol)
      }
    
      if('GOCC' %in% pathway_database)
      {
        gse_GOCC[[i-2]] <- gseGO(geneList=gene_list, ont='CC',
                                 OrgDb=org.Hs.eg.db, keyType='ENTREZID',
                                 pvalueCutoff = 1, verbose = FALSE,
                                 pAdjustMethod = "fdr")
        gene_symbol <- c()
        for(j in 1:nrow(gse_GOCC[[i-2]]))
        {
          gene_symbol <- c(gene_symbol,entrez_to_symbol(gse_GOCC[[i-2]][j,11], geneid))
        }
        gse_GOCC[[i-2]] <- cbind.data.frame(gse_GOCC[[i-2]], GeneID2=gene_symbol)
      }
    }
    for(i in pathway_database)
    {
      res_gse <- get(paste0("gse_",i))
      saveRDS(res_gse,file=paste0(output,'gse',i,'.rds'))
    }
  }
}



