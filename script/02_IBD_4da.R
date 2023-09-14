
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(GEOquery)
library(Seurat)
library(cowplot)


####################################
## data loading and preprocessing ##
####################################

input <- '/.../Data/IBD/'
GEOs <- c("GSE12251","GSE73661","GSE23597","GSE16879")

#Load gset and annot files to the environment
for(i in GEOs){
  load(paste0(input,i,".gset"))
  assign(paste0("gset",i), gset)
  load(paste0(input,i,".annot"))
  assign(paste0("annot",i),annot) }
#ADJUST GENE NAMES AND REMOVE DUPLICATED GENES
for(i in GEOs){
  annot <- get(paste0("annot",i))
  annot <- annot[annot$`Gene symbol` != "",]
  annot <- annot[grep("///",annot$`Gene symbol`,invert = T),]
  data <- exprs(get(paste0("gset",i)))[as.character(annot$ID),]
  rownames(data) <- annot$`Gene symbol`
  data <- data[!duplicated(rownames(data)),]
  data <- data[order(rownames(data)),]
  assign(paste0("data",i),data) }

#tansform all data to LOG scale
#mypar(1,4)
for (i in GEOs){
  data <- get(paste0("data",i))
  if (max(data) > 50){ data <- log(data+1,2) }
  filter <- apply(data,1,function(x){ sum(x>6) > ncol(data)*0.2 })
  data <- data[filter,]
  #  shist(data,xlim=c(5,17),main=i,col="red")
  assign(paste0("data",i),data) }
#Merge datasets
genes_in_common <- c()
for(i in GEOs){
  if(length(genes_in_common)==0){ genes_in_common <- rownames(get(paste0("data",i)))
  } else { genes_in_common <- genes_in_common[genes_in_common %in% rownames(get(paste0("data",i)))] } }
data <- matrix()
for(i in GEOs){
  if(dim(data)[1] == 1){ data <- get(paste0("data",i))[genes_in_common,]
  } else { data <- cbind(data,get(paste0("data",i))[genes_in_common,]) }}

phenotData <- read.csv(paste0(input,'merged_phenodata.csv'),row.names=1)
use_groups <- (phenotData$'week..w..ch1' %in% c("W0")) 
Group <- factor(phenotData$resp[use_groups])
phenotData <- phenotData[use_groups,c(1,11)]
colnames(phenotData)[2] <- 'Group'
data <- data[,rownames(phenotData)]

anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(GEOs))
{
  temp_anno <- phenotData[phenotData$GEO %in% GEOs[i],]
  temp_da <- data[,rownames(temp_anno)]
  temp_da_s <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(temp_da_s))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),anno_all)
colnames(anno_all)[2] <- 'DataSet'
output <- '/.../Results/IBD/hg/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg.csv'))



##############
## multiCCA ##
##############
source('/.../Code/JECCA.R')
input <- '/.../Results/IBD/hg/'
da <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg.csv'))
out <- JECCA(da,meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding.rds'))


####################
## seurat cluster ##
####################

input <- '/.../Results/IBD/hg/'
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
anno_all <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_embedding_all <- readRDS(paste0(input,'embedding.rds'))
da_all <- da_all[,rownames(anno_all)]
da_embedding_all <- da_embedding_all[rownames(anno_all),1:10]
da_seurat <- CreateSeuratObject(counts=da_all)
hvgs <- rownames(da_seurat)
da_seurat <- ScaleData(da_seurat)
da_seurat <-  RunPCA(da_seurat, features = hvgs, npcs = 30)
da_seurat[['cca']] <- CreateDimReducObject(embeddings = da_embedding_all[colnames(da_seurat),],
                                           key='CCA_', assay='RNA')
da_seurat@meta.data <- cbind.data.frame(da_seurat@meta.data, 
                                        anno_all[colnames(da_seurat),])

da_seurat <- da_seurat %>% 
  RunUMAP(reduction = "cca", dims = 1:ncol(da_embedding_all)) %>% 
  FindNeighbors(reduction = "cca", dims = 1:ncol(da_embedding_all),graph.name="rna.snn.h")

for(i in 1:20)
{
  da_seurat <- FindClusters(da_seurat, resolution = i/10,graph.name="rna.snn.h")
}
saveRDS(da_seurat, file=paste0(input,'da_seurat.rds'))

table(da_seurat$rna.snn.h_res.0.5, da_seurat$Group)
#NR  R
#0 18 36
#1 31 17
chisq.test(da_seurat$rna.snn.h_res.0.5, da_seurat$Group)
#X-squared = 8.7292, df = 1, p-value = 0.003131

###########################################
## DEG analysis between the two clusters ##
###########################################
library(limma)
anno_all <- cbind(anno_all[colnames(da_seurat),],cid=da_seurat$rna.snn.h_res.0.5)
da_all <- da_all[,rownames(anno_all)]
design <- model.matrix(~0+cid+DataSet,anno_all)
fit <- lmFit(da_all, design)
group_id1 <- 'cid1'
group_id2 <- 'cid0'
contrast_str <- paste(group_id1, group_id2, sep='-') %>% paste(collapse = ',')
prestr="makeContrasts("
poststr=",levels=design)"
commandstr=paste(prestr,contrast_str,poststr,sep="")
contrast.matrix <- eval(parse(text=commandstr))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

comp <- colnames(contrast.matrix)
deg_summary <- deg_summary2 <- out <- out_deg <- out_deg2 <- c()
for(i in 1:length(comp))
{
  temp_res <- topTable(fit2, coef=comp[i], adjust="BH", n=nrow(da_all))
  res_all <- cbind(gene=rownames(temp_res), temp_res[,-c(3,6)],comparison=comp[i])
  colnames(res_all)<-c("gene","logfc","ave_expr", "p_value","adj_p_value","comparison")
  res_deg <- filter(res_all, abs(logfc) > 1 & adj_p_value < 0.05)
  temp_counts <- c(nrow(res_deg[res_deg$logfc > 0,]),
                   nrow(res_deg[res_deg$logfc < 0,]),nrow(res_deg))
  res_deg2 <- filter(res_all, adj_p_value < 0.05)
  temp_counts2 <- c(nrow(res_deg2[res_deg2$logfc > 0,]),
                    nrow(res_deg2[res_deg2$logfc < 0,]), nrow(res_deg2))
  deg_summary <- rbind(deg_summary, temp_counts)
  deg_summary2 <- rbind(deg_summary2, temp_counts2)
  out <- rbind(out,res_all)
  out_deg <- rbind(out_deg, res_deg)
  out_deg2 <- rbind(out_deg2, res_deg2)
}
rownames(deg_summary) <- rownames(deg_summary2) <- comp
colnames(deg_summary) <- colnames(deg_summary2) <- c('Up', 'Down', 'Total')
write.csv(out, file=paste0(input,'deg_res_all.csv'))
write.csv(out_deg, file=paste0(input,'deg_fc2_fdr0.05.csv'))
write.csv(deg_summary, file=paste0(input,'deg_summary_fc2_fdr0.05.csv'))
write.csv(out_deg2, file=paste0(input,'deg_fdr0.05.csv'))
write.csv(deg_summary2, file=paste0(input,'deg_summary_fdr0.05.csv'))


###########################
## GSEA pathway analysis ##
###########################
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(AWFisher)
library(dorothea)
library(reshape2)
library(fgsea)
library(org.Hs.eg.db)

entrez_to_symbol <- function(gene, anno)
{
  gene <- gene
  anno <- as.data.frame(anno)
  entrz <- unlist(strsplit(gene, '/'))
  symbol <- anno[match(entrz, anno$entrez),]$gene
  return(paste(symbol,collapse='/'))
}

input <- '/.../Results/IBD/hg/'
da_embedding_all <- readRDS(paste0(input,'embedding.rds'))
da_embedding_all <- da_embedding_all[,1:10]
da_exprs <- read.csv(paste0(input,'exprs_scale_hg.csv'),row.names = 1)
da_exprs <- scale(da_exprs)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names=1)
GEOs <- unique(meta$DataSet)


## t-test between R/NR and CCA1-10 ##

da <- da_embedding_all[rownames(meta),]
da <- cbind.data.frame(da, Group=meta$Group)
da$Group <- factor(da$Group, levels = c('R','NR'))
t_summary <- c()
for(i in 1:(ncol(da)-1))
{
  temp <- t.test(da[,i] ~ da$Group)
  temp2 <- c(temp$statistic, temp$p.value, temp$estimate)
  t_summary <- rbind.data.frame(t_summary, temp2)
}
rownames(t_summary) <- paste0('JECCA_',1:nrow(t_summary))
colnames(t_summary) <- c('tstat', 'pvalue', 'R_mean', 'NR_mean')
write.csv(t_summary, file=paste0(input,'t_summary.csv'))


## generate expression profile for each JECCA per dataset ##
subtype_exprs <- c()
jecca_names <- paste0('JECCA',1:ncol(da_embedding_all))
for(i in 1:length(GEOs))
{
  temp_exprs <- as.matrix(da_exprs[,rownames(meta[meta$DataSet == GEOs[i],])]) %*%
    as.matrix(da_embedding_all[rownames(meta[meta$DataSet == GEOs[i],]),])
  colnames(temp_exprs) <- paste0(GEOs[i],'_',jecca_names)
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
write.csv(subtype_exprs,file=paste0(input,'subtype_exprs.csv'),row.names=F)

## get reactomPA database
ra_pathways <- reactomePathways(as.character(subtype_exprs$entrez))

## GSEA based on GO, ReactomPA ##
gse <- gse2 <- NULL
for(i in 3:ncol(subtype_exprs))
{
  gene_list <- subtype_exprs[,i]
  # name the vector
  names(gene_list) <- subtype_exprs$entrez
  # omit any NA values 
  gene_list<-na.omit(gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  gse[[i-2]] <- fgsea(ra_pathways, gene_list, 
                       minSize=10, maxSize=500)
  gse2[[i-2]] <- gseGO(geneList=gene_list, 
                       ont='CC',
                       OrgDb=org.Hs.eg.db,
                       keyType='ENTREZID',
                       pvalueCutoff = 1, 
                       verbose = FALSE,
                       pAdjustMethod = "fdr")
  gene_symbol <- gene_symbol2 <- c()
  for(j in 1:nrow(gse[[i-2]]))
  {
    gse[[i-2]]$leadingEdge[j] <- paste(unlist(gse[[i-2]]$leadingEdge[j]),
                                        collapse='/')
    gene_symbol <- c(gene_symbol,entrez_to_symbol(unlist(gse[[i-2]][j,8]), geneid))
  }
  for(j in 1:nrow(gse2[[i-2]]))
  {
    gene_symbol2 <- c(gene_symbol2,entrez_to_symbol(gse2[[i-2]][j,11], geneid))
  }
  
  gse[[i-2]] <- cbind.data.frame(gse[[i-2]], GeneID2=gene_symbol)
  gse2[[i-2]] <- cbind.data.frame(gse2[[i-2]], GeneID2=gene_symbol2)
  if((i-2)%%10 == 0) cat('#')
}
saveRDS(gse, paste0(input,'gseReactomPA.rds'))
saveRDS(gse2, paste0(input,'gseGO.rds'))

