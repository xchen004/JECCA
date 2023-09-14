
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(irlba)
library(GEOquery)
library(PMA)
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

## add GSE73661 VDZ patients ##
phenotData2 <- read.csv(paste0(input,'GSE73661_pData.csv'),row.names = 1)
phenotData2 <- phenotData2[grep('vdz',phenotData2$induction.therapy_maintenance.therapy.ch1),]
phenotData2 <- phenotData2[,c(1,36:38,40)]
temp_pos <- regexpr('_',phenotData2[,1])
phenotData2$resp <- substr(phenotData2[,1],temp_pos+1,999)
temp_pos <- regexpr('_',phenotData2$resp)
phenotData2$resp <- substr(phenotData2$resp,1,temp_pos-1)
phenotData2$resp[grep(' R ',phenotData2$resp)] <- 'R'
phenotData2$resp[grep(' NR ',phenotData2$resp)] <- 'NR'
temp_w0 <- phenotData2[phenotData2$week..w..ch1 %in% 'W0',-c(1,5,6)]
temp_w0 <- cbind.data.frame(GEO_id=rownames(temp_w0), temp_w0)
temp_w6 <- phenotData2[phenotData2$week..w..ch1 %in% 'W6',-c(1,5)]
temp_w12 <- phenotData2[phenotData2$week..w..ch1 %in% 'W12',-c(1,5)]
temp_w52 <- phenotData2[phenotData2$week..w..ch1 %in% 'W52',-c(1,5)]
colnames(temp_w0)[-1] <- c('Trt', 'mayo_score_BL','sid')
colnames(temp_w6) <- c('Trt', 'mayo_score_w6','sid','resp_w6')
colnames(temp_w12) <- c('Trt', 'mayo_score_w12','sid','resp_w12')
colnames(temp_w52) <- c('Trt', 'mayo_score_w52','sid','resp_w52')
phenotData2 <- merge(temp_w0, temp_w6, by=c('sid','Trt'),all=T)
phenotData2 <- merge(phenotData2, temp_w12, by=c('sid','Trt'),all=T)
phenotData2 <- merge(phenotData2, temp_w52, by=c('sid','Trt'),all=T)
phenotData2[is.na(phenotData2)] <- ''
phenotData2$GEO <- 'GSE73661_VDZ'
phenotData2$Group <- factor(apply(phenotData2[,grep('resp',colnames(phenotData2))],1,paste,collapse=''))
levels(phenotData2$Group) <- c('NR','NR','R','R','R')
rownames(phenotData2) <- phenotData2$GEO_id
phenotData2 <- phenotData2[,11:12]
phenotData <- rbind.data.frame(phenotData, phenotData2)
data <- data[,rownames(phenotData)]

## add GSE109142 ##
library(edgeR)
library(sva)
GSE109142_meta <- pData(getGEO("GSE109142",GSEMatrix =TRUE)[[1]])
GSE109142_meta <- GSE109142_meta[GSE109142_meta$`diagnosis:ch1` == "Ulcerative Colitis",]
dataGSE109142 <-  read.csv(paste0(input,'rawdataGSE109142.csv'),row.names = 1)
#Merge data and remove 0 counts
data2 <- as.matrix(dataGSE109142)
sel <- apply(data2,1,function(x) sum(x > 5) > ncol(data2)*0.1 )
data2 <- data2[sel,]
#Merge metadata
phenotData3 <- GSE109142_meta
Group <- phenotData3$`week 4 remission:ch1`
## remove sex effect ##
modCombat <- model.matrix(~1,data=phenotData3)
combat_edata <- ComBat(dat=data2,batch=phenotData3$`Sex:ch1`,mod=modCombat)
combat_edata[combat_edata<0] <- 0
y <- DGEList(counts=combat_edata, modCombat,remove.zeros = T)
y <- calcNormFactors(y,method = "TMM")
y <- estimateCommonDisp(y)
da2 <- cpm(y, log = T, prior.count = 1, normalized.lib.sizes = T)
phenotData3 <- cbind.data.frame(GEO_id=rownames(phenotData3),
                                GEO='GSE109142',
                                Group=factor(phenotData3$`week 4 remission:ch1`))
levels(phenotData3$Group) <- c('NR','R')
rownames(phenotData3) <- phenotData3$GEO_id
phenotData3 <- phenotData3[,-1]
da2 <- da2[,rownames(phenotData3)]
common_genes <- intersect(rownames(data),rownames(da2))
data <- cbind(data[common_genes,], da2[common_genes,])
phenotData <- rbind.data.frame(phenotData, phenotData3)

GEOs <- c(GEOs, 'GSE73661_VDZ', 'GSE109142')
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
output <- '/.../Results/IBD/hg_6da/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg.csv'))



##############
## multiCCA ##
##############
source('/.../Code/JECCA.R')
input <- '/.../Results/IBD/hg_6da/'
da <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg.csv'))
out <- JECCA(da,meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding.rds'))


#############################
## classification analysis ##
#############################

library(precrec)
rownames(meta) <- meta$sid
ds <- names(table(meta$DataSet))
test_geo <- c('GSE109142', 'GSE73661_VDZ')
train_geo <- ds[2:5]
da <- readRDS(paste0(input,'embedding.rds'))
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet %in% train_geo & !is.na(meta$Group),]
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
  temp_stat <- c(auc(res_roc[[i]])[1,4])
  res_stat <- c(res_stat, temp_stat)
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred.csv'))


###########
## ComBat##
###########

library(sva)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
dir.create(output, recursive = T)
modCombat <- model.matrix(~1,data=meta)
da_combat <- ComBat(dat=da_all,batch=meta$DataSet,mod=modCombat)

da_seurat <- CreateSeuratObject(counts=da_combat)
hvgs <- rownames(da_seurat)
da_seurat <- ScaleData(da_seurat)
da_seurat <-  RunPCA(da_seurat, features = hvgs, npcs = 30)
da_seurat@meta.data <- cbind.data.frame(da_seurat@meta.data, 
                                        meta[match(colnames(da_seurat), rownames(meta)),])

da_seurat <- da_seurat %>% RunUMAP(reduction = "pca",dims=1:30) %>% 
  FindNeighbors(reduction = "pca",graph.name="rna.snn.h",dims=1:30)
da_embedding <- Embeddings(da_seurat, reduction = 'pca')
saveRDS(da_embedding, file=paste0(output,'embedding_combat.rds'))


########################################
## classification analysis for ComBat ##
########################################

ds <- names(table(meta$DataSet))
test_geo <- c('GSE109142', 'GSE73661_VDZ')
train_geo <- ds[2:5]
da <- readRDS(paste0(input,'embedding_combat.rds'))
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet %in% train_geo & !is.na(meta$Group),]
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
  temp_stat <- c(auc(res_roc[[i]])[1,4])
  res_stat <- c(res_stat, temp_stat)
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_combat.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_combat.csv'))

