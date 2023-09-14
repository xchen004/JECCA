
#########
## DSS ##
#########

####################################
## data loading and preprocessing ##
####################################

input <- '/.../Results/IBD/hg/'
input2 <- '/.../Data/Mouse/'
data <- read.csv(paste0(input,'exprs_hg.csv'),row.names = 1)
phenotData <- read.csv(paste0(input,'anno_hg.csv'),row.names=1)
da_mouse <- read.csv(paste0(input2,'GSE131032_log2_counts_per_million.csv'),row.names = 1)
temp <- read.csv(paste0(input2,'metadata_DSS.csv'))
phenotData2 <- cbind(DataSet='GSE131032',Group = temp$Group)
rownames(phenotData2) <- temp$SampleID
phenotData <- rbind(phenotData, phenotData2)
da_mouse <- da_mouse[,rownames(phenotData2)]
homolog <- read.table(paste0(input2,'homologs.txt'),header=T,sep=',',row.names = 1)

data <- cbind.data.frame(hg_gene=rownames(data),data)
da_mouse <- cbind.data.frame(ms_gene=rownames(da_mouse),da_mouse)
colnames(homolog) <- c('ms_gene','hg_gene')
data <- merge(data,homolog,by='hg_gene')
data <- merge(data, da_mouse, by='ms_gene')
rownames(data) <- paste0(data$hg_gene,'_',data$ms_gene)
data <- data[,-c(1:2)]
GEOs <- unique(phenotData$DataSet)

da_s <- anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(GEOs))
{
  temp_anno <- phenotData[phenotData$DataSet %in% GEOs[i],]
  temp_da <- data[,rownames(temp_anno)]
  da_s[[i]] <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(da_s[[i]]))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),anno_all)
output <- '/.../Results/Cross_species/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg_DSS.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg_DSS.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg_DSS.csv'))


##############
## multiCCA ##
##############
source('/.../Code/JECCA.R')
input <- '/.../Results/Cross_species/'
da <- read.csv(paste0(input,'exprs_hg_DSS.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg_DSS.csv'))
out <- JECCA(da,meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding_DSS.rds'))


## t-test between R/NR and CCA1-10 ##
input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_DSS.rds'))
da_embedding_all <- da_embedding_all[,1:10]
meta <- read.csv(paste0(input,'anno_hg_DSS.csv'),row.names=1)
meta <- meta[-grep('day',meta$Group),]
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
write.csv(t_summary, file=paste0(input,'t_summary_DSS.csv'))


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

input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_DSS.rds'))
da_embedding_all <- da_embedding_all[,1:10]
da_exprs <- read.csv(paste0(input,'exprs_scale_hg_DSS.csv'),row.names = 1)
da_exprs <- scale(da_exprs)
meta <- read.csv(paste0(input,'anno_hg_DSS.csv'),row.names=1)
da_exprs <- da_exprs[,rownames(meta)]
GEOs <- unique(meta$DataSet)
da_exprs <- cbind.data.frame(gene=matrix(unlist(strsplit(rownames(da_exprs),'_')),nrow=2)[1,], 
                             da_exprs)
da_exprs <- da_exprs[!duplicated(da_exprs$gene),]
rownames(da_exprs) <- da_exprs$gene
da_exprs <- da_exprs[,-1]
da_embedding_all <- da_embedding_all[rownames(meta),]


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
write.csv(subtype_exprs,file=paste0(input,'subtype_exprs_DSS.csv'),row.names=F)

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
saveRDS(gse, paste0(input,'gseReactomPA_DSS.rds'))
saveRDS(gse2, paste0(input,'gseGO_DSS.rds'))


#############################################################
## classification to get probability of R/NR for DSS mouse ##
#############################################################

input <- '/.../Results/Cross_species/'
meta <- read.csv(paste0(input,'anno_hg_DSS.csv'),row.names=1)
da <- readRDS(paste0(input,'embedding_DSS.rds'))
da <- da[rownames(meta),]
meta_train <- meta[meta$DataSet != 'GSE131032' & !is.na(meta$Group),]
meta_test <- meta[meta$DataSet == 'GSE131032' & !is.na(meta$Group),]
da_train <- da[rownames(meta_train),]
da_test <- da[rownames(meta_test),]
da_train2 <- cbind.data.frame(resp=factor(meta_train$Group), da_train)
da_test2 <- cbind.data.frame(resp=factor(meta_test$Group), da_test)
set.seed(123)
train_ranger <- ranger::ranger(resp ~ ., data = da_train2, probability = TRUE)
pred_ranger <- predict(train_ranger, da_test2)
pred_prob <- cbind.data.frame(Group=meta_test$Group, 
                              round(pred_ranger$prediction,3))
write.csv(pred_prob, file=paste0(input,'DSS_mouse_R_NR_prediction.csv'),row.names=F)

#####################################################
## classification to get day status for human data ##
#####################################################

set.seed(123)
train_ranger <- ranger::ranger(resp ~ ., data = da_test2)
pred_ranger <- predict(train_ranger, da_train2)
pred <- cbind.data.frame(Group=meta_train$Group, prediction=pred_ranger$prediction,
                         prediction2=factor(pred_ranger$prediction))
levels(pred$prediction2) <- paste0('day',c('00','02-04','02-04','06-08','06-08',
                                           '06-08','10-12','10-12','14'))
write.csv(pred, file=paste0(input,'DSS_human_day0_day14_prediction.csv'),row.names=F)


###############################
## classification to get AUC ##
###############################

library(precrec)
ds <- unique(meta$DataSet)
test_geo <- ds[-5]
da <- readRDS(paste0(input,'embedding_DSS.rds'))
meta <- meta[meta$DataSet %in% ds[-5],]
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:10])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:10])
  
  set.seed(123)
  train_ranger <- ranger::ranger(Group ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  nr_prob <- pred_ranger$prediction[,1]
  pred_Group <- rep('NR',nrow(da_test2))
  pred_Group[nr_prob < 0.5] <- 'R'
  temp <- cbind.data.frame(Group=da_test2$Group, pred_Group=pred_Group, 
                           nr_prob=nr_prob, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  res <- rbind.data.frame(res, temp)
  
  res_roc[[i]] <- evalmod(scores = 1-temp$nr_prob, labels = temp$Group)
  df <- cbind.data.frame(Reference=factor(temp$Group), 
                         Prediction=factor(temp$pred_Group))
  res_stat <- c(res_stat, auc(res_roc[[i]])[1,4])
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_DSS.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_DSS.csv'))


#########
## TNF ##
#########

####################################
## data loading and preprocessing ##
####################################
library(edgeR)
library(biomaRt)
input <- '/.../Results/IBD/hg/'
input2 <- '/.../Data/Mouse/'
data <- read.csv(paste0(input,'exprs_hg.csv'),row.names = 1)
phenotData <- read.csv(paste0(input,'anno_hg.csv'),row.names=1)
da_tct <- read.table(paste0(input2,'GSE164339_TCT_Raw_gene_counts_matrix.txt'),
                     row.names = 1, header=T)
da_tnf <- read.table(paste0(input2,'GSE164339_TNF_Raw_gene_counts_matrix.txt'),
                     row.names = 1, header=T)
da_mouse <- cbind(da_tct, da_tnf[rownames(da_tct),])
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
geneid <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                filters = 'ensembl_gene_id', values = rownames(da_mouse), mart = ensembl)
geneid2 <- geneid[!duplicated(geneid[,2]),]
da_mouse <- da_mouse[geneid2[,1],]
rownames(da_mouse) <- geneid2[,2]
phenotData2 <- read.csv(paste0(input2,'metadata_TNF_TCT.csv'),row.names = 1)
phenotData2$DataSet <- 'GSE164339'
phenotData2 <- phenotData2[phenotData2$Treatment %in% 'Vehicle' & phenotData2$strain %in% 'TNF',]
phenotData2 <- phenotData2[,c(4,2)]
colnames(phenotData2)[2] <- 'Group'
phenotData <- rbind(phenotData, phenotData2)
da_mouse <- da_mouse[,rownames(phenotData2)]
y <- DGEList(counts=da_mouse, remove.zeros = T)
y <- calcNormFactors(y,method = "TMM")
y <- estimateCommonDisp(y)
da_mouse <- cpm(y, log = T, prior.count = 1, normalized.lib.sizes = T)

homolog <- read.table(paste0(input2,'homologs.txt'),header=T,sep=',',row.names = 1)
colnames(homolog) <- c('ms_gene','hg_gene')
data <- cbind.data.frame(hg_gene=rownames(data),data)
da_mouse <- cbind.data.frame(ms_gene=rownames(da_mouse),da_mouse)
data <- merge(data,homolog,by='hg_gene')
data <- merge(data, da_mouse, by='ms_gene')
rownames(data) <- paste0(data$hg_gene,'_',data$ms_gene)
data <- data[,-c(1:2)]
GEOs <- unique(phenotData$DataSet)

da_s <- anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(GEOs))
{
  temp_anno <- phenotData[phenotData$DataSet %in% GEOs[i],]
  temp_da <- data[,rownames(temp_anno)]
  da_s[[i]] <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(da_s[[i]]))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),anno_all)
output <- '/.../Results/Cross_species/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg_TNF.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg_TNF.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg_TNF.csv'))


###############
## multiCCA  ##
###############
source('/.../Code/JECCA.R')
input <- '/.../Results/Cross_species/'
da <- read.csv(paste0(input,'exprs_hg_TNF.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg_TNF.csv'))
out <- JECCA(da,meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding_TNF.rds'))


## t-test between R/NR and CCA1-10 ##
input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_TNF.rds'))
da_embedding_all <- da_embedding_all[,1:10]
meta <- read.csv(paste0(input,'anno_hg_TNF.csv'),row.names=1)
meta <- meta[-grep('flame',meta$Group),]
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
write.csv(t_summary, file=paste0(input,'t_summary_TNF.csv'))



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

input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_TNF.rds'))
da_embedding_all <- da_embedding_all[,1:10]
da_exprs <- read.csv(paste0(input,'exprs_scale_hg_TNF.csv'),row.names = 1)
da_exprs <- scale(da_exprs)
meta <- read.csv(paste0(input,'anno_hg_TNF.csv'),row.names=1)
da_exprs <- da_exprs[,rownames(meta)]
GEOs <- unique(meta$DataSet)
da_exprs <- cbind.data.frame(gene=matrix(unlist(strsplit(rownames(da_exprs),'_')),nrow=2)[1,], 
                             da_exprs)
da_exprs <- da_exprs[!duplicated(da_exprs$gene),]
rownames(da_exprs) <- da_exprs$gene
da_exprs <- da_exprs[,-1]
da_embedding_all <- da_embedding_all[rownames(meta),]


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
write.csv(subtype_exprs,file=paste0(input,'subtype_exprs_TNF.csv'),row.names=F)

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
saveRDS(gse, paste0(input,'gseReactomPA_TNF.rds'))
saveRDS(gse2, paste0(input,'gseGO_TNF.rds'))


#############################################################
## classification to get probability of R/NR for TNF mouse ##
#############################################################
input <- '/.../Results/Cross_species/'
meta <- read.csv(paste0(input,'anno_hg_TNF.csv'),row.names=1)
da <- readRDS(paste0(input,'embedding_TNF.rds'))
da <- da[rownames(meta),]
meta_train <- meta[meta$DataSet != 'GSE164339' & !is.na(meta$Group),]
meta_test <- meta[meta$DataSet == 'GSE164339' & !is.na(meta$Group),]
da_train <- da[rownames(meta_train),]
da_test <- da[rownames(meta_test),]
da_train2 <- cbind.data.frame(resp=factor(meta_train$Group), da_train)
da_test2 <- cbind.data.frame(resp=factor(meta_test$Group), da_test)
set.seed(123)
train_ranger <- ranger::ranger(resp ~ ., data = da_train2, probability = TRUE)
pred_ranger <- predict(train_ranger, da_test2)
pred_prob <- cbind.data.frame(Group=meta_test$Group, 
                              round(pred_ranger$prediction,3))
write.csv(pred_prob, file=paste0(input,'TNF_mouse_R_NR_prediction.csv'),row.names=F)


###############################
## classification to get AUC ##
###############################

library(precrec)
ds <- unique(meta$DataSet)
test_geo <- ds[-5]
da <- readRDS(paste0(input,'embedding_TNF.rds'))
meta <- meta[meta$DataSet %in% ds[-5],]
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:10])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:10])
  
  set.seed(123)
  train_ranger <- ranger::ranger(Group ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  nr_prob <- pred_ranger$prediction[,1]
  pred_Group <- rep('NR',nrow(da_test2))
  pred_Group[nr_prob < 0.5] <- 'R'
  temp <- cbind.data.frame(Group=da_test2$Group, pred_Group=pred_Group, 
                           nr_prob=nr_prob, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  res <- rbind.data.frame(res, temp)
  
  res_roc[[i]] <- evalmod(scores = 1-temp$nr_prob, labels = temp$Group)
  df <- cbind.data.frame(Reference=factor(temp$Group), 
                         Prediction=factor(temp$pred_Group))
  res_stat <- c(res_stat, auc(res_roc[[i]])[1,4])
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_TNF.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_TNF.csv'))


#########
## TCT ##
#########

####################################
## data loading and preprocessing ##
####################################
library(edgeR)
library(biomaRt)
input <- '/.../Results/IBD/hg/'
input2 <- '/.../Data/Mouse/'
data <- read.csv(paste0(input,'exprs_hg.csv'),row.names = 1)
phenotData <- read.csv(paste0(input,'anno_hg.csv'),row.names=1)
da_tct <- read.table(paste0(input2,'GSE164339_TCT_Raw_gene_counts_matrix.txt'),
                     row.names = 1, header=T)
da_tnf <- read.table(paste0(input2,'GSE164339_TNF_Raw_gene_counts_matrix.txt'),
                     row.names = 1, header=T)
da_mouse <- cbind(da_tct, da_tnf[rownames(da_tct),])
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
geneid <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                filters = 'ensembl_gene_id', values = rownames(da_mouse), mart = ensembl)
geneid2 <- geneid[!duplicated(geneid[,2]),]
da_mouse <- da_mouse[geneid2[,1],]
rownames(da_mouse) <- geneid2[,2]
phenotData2 <- read.csv(paste0(input2,'metadata_TNF_TCT.csv'),row.names = 1)
phenotData2$DataSet <- 'GSE164339'
phenotData2 <- phenotData2[phenotData2$Treatment %in% 'Vehicle' & phenotData2$strain %in% 'TCT',]
phenotData2 <- phenotData2[,c(4,2)]
colnames(phenotData2)[2] <- 'Group'
phenotData <- rbind(phenotData, phenotData2)
da_mouse <- da_mouse[,rownames(phenotData2)]
y <- DGEList(counts=da_mouse, remove.zeros = T)
y <- calcNormFactors(y,method = "TMM")
y <- estimateCommonDisp(y)
da_mouse <- cpm(y, log = T, prior.count = 1, normalized.lib.sizes = T)

homolog <- read.table(paste0(input2,'homologs.txt'),header=T,sep=',',row.names = 1)
colnames(homolog) <- c('ms_gene','hg_gene')
data <- cbind.data.frame(hg_gene=rownames(data),data)
da_mouse <- cbind.data.frame(ms_gene=rownames(da_mouse),da_mouse)
data <- merge(data,homolog,by='hg_gene')
data <- merge(data, da_mouse, by='ms_gene')
rownames(data) <- paste0(data$hg_gene,'_',data$ms_gene)
data <- data[,-c(1:2)]
GEOs <- unique(phenotData$DataSet)

da_s <- anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(GEOs))
{
  temp_anno <- phenotData[phenotData$DataSet %in% GEOs[i],]
  temp_da <- data[,rownames(temp_anno)]
  da_s[[i]] <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(da_s[[i]]))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),anno_all)
output <- '/.../Results/Cross_species/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg_TCT.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg_TCT.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg_TCT.csv'))


###############
## multiCCA  ##
###############
source('/.../Code/JECCA.R')
input <- '/.../Results/Cross_species/'
da <- read.csv(paste0(input,'exprs_hg_TCT.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg_TCT.csv'))
out <- JECCA(da,meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding_TCT.rds'))


## t-test between R/NR and CCA1-10 ##
input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_TCT.rds'))
da_embedding_all <- da_embedding_all[,1:10]
meta <- read.csv(paste0(input,'anno_hg_TCT.csv'),row.names=1)
meta <- meta[-grep('flame',meta$Group),]
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
write.csv(t_summary, file=paste0(input,'t_summary_TCT.csv'))



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

input <- '/.../Results/Cross_species/'
da_embedding_all <- readRDS(paste0(input,'embedding_TCT.rds'))
da_embedding_all <- da_embedding_all[,1:10]
da_exprs <- read.csv(paste0(input,'exprs_scale_hg_TCT.csv'),row.names = 1)
da_exprs <- scale(da_exprs)
meta <- read.csv(paste0(input,'anno_hg_TCT.csv'),row.names=1)
da_exprs <- da_exprs[,rownames(meta)]
GEOs <- unique(meta$DataSet)
da_exprs <- cbind.data.frame(gene=matrix(unlist(strsplit(rownames(da_exprs),'_')),nrow=2)[1,], 
                             da_exprs)
da_exprs <- da_exprs[!duplicated(da_exprs$gene),]
rownames(da_exprs) <- da_exprs$gene
da_exprs <- da_exprs[,-1]
da_embedding_all <- da_embedding_all[rownames(meta),]


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
write.csv(subtype_exprs,file=paste0(input,'subtype_exprs_TCT.csv'),row.names=F)

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
saveRDS(gse, paste0(input,'gseReactomPA_TCT.rds'))
saveRDS(gse2, paste0(input,'gseGO_TCT.rds'))


#######################################################
## classification to get probability of R/NR for TCT ##
#######################################################
input <- '/.../Results/Cross_species/'
meta <- read.csv(paste0(input,'anno_hg_TCT.csv'),row.names=1)
da <- readRDS(paste0(input,'embedding_TCT.rds'))
da <- da[rownames(meta),]
meta_train <- meta[meta$DataSet != 'GSE164339' & !is.na(meta$Group),]
meta_test <- meta[meta$DataSet == 'GSE164339' & !is.na(meta$Group),]
da_train <- da[rownames(meta_train),]
da_test <- da[rownames(meta_test),]
da_train2 <- cbind.data.frame(resp=factor(meta_train$Group), da_train)
da_test2 <- cbind.data.frame(resp=factor(meta_test$Group), da_test)
set.seed(123)
train_ranger <- ranger::ranger(resp ~ ., data = da_train2, probability = TRUE)
pred_ranger <- predict(train_ranger, da_test2)
pred_prob <- cbind.data.frame(Group=meta_test$Group, 
                              round(pred_ranger$prediction,3))
write.csv(pred_prob, file=paste0(input,'TCT_mouse_R_NR_prediction.csv'),row.names=F)


###############################
## classification to get AUC ##
###############################

library(precrec)
ds <- unique(meta$DataSet)
test_geo <- ds[-5]
da <- readRDS(paste0(input,'embedding_TCT.rds'))
meta <- meta[meta$DataSet %in% ds[-5],]
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:10])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:10])
  
  set.seed(123)
  train_ranger <- ranger::ranger(Group ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  nr_prob <- pred_ranger$prediction[,1]
  pred_Group <- rep('NR',nrow(da_test2))
  pred_Group[nr_prob < 0.5] <- 'R'
  temp <- cbind.data.frame(Group=da_test2$Group, pred_Group=pred_Group, 
                           nr_prob=nr_prob, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  res <- rbind.data.frame(res, temp)
  
  res_roc[[i]] <- evalmod(scores = 1-temp$nr_prob, labels = temp$Group)
  df <- cbind.data.frame(Reference=factor(temp$Group), 
                         Prediction=factor(temp$pred_Group))
  res_stat <- c(res_stat, auc(res_roc[[i]])[1,4])
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_TCT.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_TCT.csv'))



## pathway comparison of three mouse models ##
input <- '/.../Results/Cross_species/'
gse_dss <- readRDS(paste0(input,'gseReactomPA_DSS.rds'))
gse_tnf <- readRDS(paste0(input,'gseReactomPA_TNF.rds'))
gse_tct <- readRDS(paste0(input,'gseReactomPA_TCT.rds'))
subtype_exprs <- read.csv(paste0(input,'subtype_exprs_DSS.csv'))
pathway_id_dss <- gse_dss[[1]][,1]
gse_dss_p <- gse_dss[[1]]$pval
gse_dss_NES <- gse_dss[[1]]$NES
pathway_id_tnf <- gse_tnf[[1]][,1]
gse_tnf_p <- gse_tnf[[1]]$pval
gse_tnf_NES <- gse_tnf[[1]]$NES
pathway_id_tct <- gse_tct[[1]][,1]
gse_tct_p <- gse_tct[[1]]$pval
gse_tct_NES <- gse_tct[[1]]$NES
for(i in 2:length(gse_dss))
{
  gse_dss_p <- cbind.data.frame(gse_dss_p, gse_dss[[i]][match(pathway_id_dss,gse_dss[[i]][,1]),]$padj)
  gse_dss_NES <- cbind.data.frame(gse_dss_NES, gse_dss[[i]][match(pathway_id_dss,gse_dss[[i]][,1]),]$NES)
  gse_tnf_p <- cbind.data.frame(gse_tnf_p, gse_tnf[[i]][match(pathway_id_tnf,gse_tnf[[i]][,1]),]$padj)
  gse_tnf_NES <- cbind.data.frame(gse_tnf_NES, gse_tnf[[i]][match(pathway_id_tnf,gse_tnf[[i]][,1]),]$NES)
  gse_tct_p <- cbind.data.frame(gse_tct_p, gse_tct[[i]][match(pathway_id_tct,gse_tct[[i]][,1]),]$padj)
  gse_tct_NES <- cbind.data.frame(gse_tct_NES, gse_tct[[i]][match(pathway_id_tct,gse_tct[[i]][,1]),]$NES)
}
colnames(gse_dss_p) <- colnames(gse_dss_NES) <- 
  colnames(gse_tnf_p) <- colnames(gse_tnf_NES) <- 
  colnames(gse_tct_p) <- colnames(gse_tct_NES) <- colnames(subtype_exprs)[-c(1:2)]
rownames(gse_dss_p) <- rownames(gse_dss_NES) <- pathway_id_dss
rownames(gse_tnf_p) <- rownames(gse_tnf_NES) <- pathway_id_tnf
rownames(gse_tct_p) <- rownames(gse_tct_NES) <- pathway_id_tct

overlap_pathways <- intersect(pathway_id_dss, intersect(pathway_id_tnf, pathway_id_tct))
gse_dss_p <- gse_dss_p[overlap_pathways,grep('GSE16879',colnames(gse_dss_p))]
gse_dss_NES <- gse_dss_NES[overlap_pathways,grep('GSE16879',colnames(gse_dss_NES))]
gse_tnf_p <- gse_tnf_p[overlap_pathways,grep('GSE16879',colnames(gse_tnf_p))]
gse_tnf_NES <- gse_tnf_NES[overlap_pathways,grep('GSE16879',colnames(gse_tnf_NES))]
gse_tct_p <- gse_tct_p[overlap_pathways,grep('GSE16879',colnames(gse_tct_p))]
gse_tct_NES <- gse_tct_NES[overlap_pathways,grep('GSE16879',colnames(gse_tct_NES))]
colnames(gse_dss_p) <- colnames(gse_dss_NES) <- sub('GSE16879','DSS',colnames(gse_dss_p))
colnames(gse_tnf_p) <- colnames(gse_tnf_NES) <- sub('GSE16879','TNF',colnames(gse_tnf_p))
colnames(gse_tct_p) <- colnames(gse_tct_NES) <- sub('GSE16879','TCT',colnames(gse_tct_p))

gse_indc <- cbind.data.frame(DSS=ifelse(rowSums(gse_dss_p < 0.05)>0, 1, 0),
                  TNF=ifelse(rowSums(gse_tnf_p < 0.05)>0, 1, 0),
                  TCT=ifelse(rowSums(gse_tct_p < 0.05)>0, 1, 0))
rownames(gse_indc) <- overlap_pathways

pathway_all <- rownames(gse_indc[rowSums(gse_indc)==3,])
pathway_all <- cbind(gse_dss_p[pathway_all,],gse_tnf_p[pathway_all,],gse_tct_p[pathway_all,])
pathway_dss_tnf <- rownames(gse_indc[rowSums(gse_indc)==2 & gse_indc$TCT==0,])
pathway_dss_tnf <- cbind(gse_dss_p[pathway_dss_tnf,],gse_tnf_p[pathway_dss_tnf,])
pathway_dss_tct <- rownames(gse_indc[rowSums(gse_indc)==2 & gse_indc$TNF==0,])
pathway_dss_tct <- cbind(gse_dss_p[pathway_dss_tct,], gse_tct_p[pathway_dss_tct,])
pathway_tnf_tct <- rownames(gse_indc[rowSums(gse_indc)==2 & gse_indc$DSS==0,])
pathway_tnf_tct <- cbind(gse_tnf_p[pathway_tnf_tct,],gse_tct_p[pathway_tnf_tct,])
pathway_dss <- gse_dss_p[rownames(gse_indc[rowSums(gse_indc)==1 & gse_indc$DSS==1,]),]
pathway_dss <- pathway_dss[order(apply(pathway_dss,1,min)),]
pathway_tnf <- gse_tnf_p[rownames(gse_indc[rowSums(gse_indc)==1 & gse_indc$TNF==1,]),]
pathway_tnf <- pathway_tnf[order(apply(pathway_tnf,1,min)),]
pathway_tct <- gse_tct_p[rownames(gse_indc[rowSums(gse_indc)==1 & gse_indc$TCT==1,]),]
pathway_tct <- pathway_tct[order(apply(pathway_tct,1,min)),]

write.csv(pathway_all,file=paste0(input,'pathway_DSS_TNF_TCT.csv'))
write.csv(pathway_dss_tnf,file=paste0(input,'pathway_DSS_TNF.csv'))
write.csv(pathway_dss_tct,file=paste0(input,'pathway_DSS_TCT.csv'))
write.csv(pathway_tnf_tct,file=paste0(input,'pathway_TNF_TCT.csv'))
write.csv(pathway_dss,file=paste0(input,'pathway_DSS.csv'))
write.csv(pathway_tnf,file=paste0(input,'pathway_TNF.csv'))
write.csv(pathway_tct,file=paste0(input,'pathway_TCT.csv'))




