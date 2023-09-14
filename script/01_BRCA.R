
####################################
## data loading and preprocessing ##
####################################

input <- '/.../Data/BRCA/'
data <- read.table(paste0(input,'KMPLOT_BRCA_XP.txt'),row.names = 1, header=T)
phenotData <- read.table(paste0(input,'KMPLOT_BRCA_SURVIVAL.txt'),
                         header=T,sep='\t',row.names=1)
## remove gse4922 which only has 1 sample, ##
## and gse12276 which doesn't have clinical information ##
ds <- names(table(phenotData$DataSet))[-c(3,9)]
anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(ds))
{
  temp_anno <- phenotData[phenotData$DataSet %in% ds[i],]
  temp_da <- data[,rownames(temp_anno)]
  temp_da_s <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(temp_da_s))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),DataSet=toupper(anno_all$DataSet),
                             ER=anno_all$ER, Grade=anno_all$Grade)
output <- '/.../Results/BRCA/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg.csv'))



##############
## multiCCA ##
##############
source('/.../Code/JECCA.R')
input <- '/.../Results/BRCA/'
da <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg.csv'))
out <- JECCA(da,meta,ncomp=30)
saveRDS(out, file=paste0(input,'embedding.rds'))


####################
## seurat cluster ##
####################
library(Seurat)
library(dplyr)
input <- '/.../Results/BRCA/'
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_embedding_all <- readRDS(paste0(input,'embedding.rds'))
da_all <- da_all[,rownames(meta)]
da_embedding_all <- da_embedding_all[rownames(meta),]
da_seurat <- CreateSeuratObject(counts=da_all)
hvgs <- rownames(da_seurat)
da_seurat <- ScaleData(da_seurat)
da_seurat <-  RunPCA(da_seurat, features = hvgs, npcs = 30)
da_seurat[['cca']] <- CreateDimReducObject(embeddings = da_embedding_all[colnames(da_seurat),],
                                           key='CCA_', assay='RNA')
da_seurat@meta.data <- cbind.data.frame(da_seurat@meta.data, 
                                        meta[colnames(da_seurat),])

da_seurat <- da_seurat %>% 
  RunUMAP(reduction = "cca", dims = 1:ncol(da_embedding_all)) %>% 
  FindNeighbors(reduction = "cca", dims = 1:ncol(da_embedding_all),graph.name="rna.snn.h")

for(i in 1:20)
{
  da_seurat <- FindClusters(da_seurat, resolution = i/10,graph.name="rna.snn.h")
}
saveRDS(da_seurat, file=paste0(input,'da_seurat.rds'))


#############################
## classification analysis ##
#############################
library(precrec)
library(Boruta)

## ER ##
da <- readRDS(paste0(input,'embedding.rds'))
ds <- names(table(meta$DataSet))
test_geo <- paste0('GSE',c(2034, 2990, 3494, 7390))
res_roc <- NULL
res <- res_stat <- c()
for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$ER),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$ER),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  set.seed(123)
  res.boruta <- Boruta(da_train, meta_train$ER)
  
  da_train2 <- cbind.data.frame(ER=factor(meta_train$ER), 
                                da_train[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  da_test2 <- cbind.data.frame(ER=factor(meta_test$ER), 
                               da_test[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  set.seed(123)
  train_ranger <- ranger::ranger(ER ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  er_prob <- pred_ranger$prediction[,1]
  pred_ER <- rep(0,nrow(da_test2))
  pred_ER[er_prob < 0.5] <- 1
  temp <- cbind.data.frame(ER=da_test2$ER, pred_ER=pred_ER, 
                           pred_p=er_prob, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  res <- rbind.data.frame(res, temp)
  
  res_roc[[i]] <- evalmod(scores = 1-temp$pred_p, labels = temp$ER)
  df <- cbind.data.frame(Reference=factor(temp$ER), Prediction=factor(temp$pred_ER))
  res_stat <- c(res_stat, auc(res_roc[[i]])[1,4])
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_ER.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_ER.csv'))


## Grade ##
library(DiagTest3Grp)
test_geo <- paste0('GSE',c(11121, 1456, 2990, 3494, 7390, 9195))
res <- res_stat <- c()
for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Grade),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Grade),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  set.seed(123)
  res.boruta <- Boruta(da_train, meta_train$Grade)
  
  da_train2 <- cbind.data.frame(Grade=factor(meta_train$Grade), 
                                da_train[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  da_test2 <- cbind.data.frame(Grade=factor(meta_test$Grade), 
                               da_test[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  set.seed(123)
  train_ranger <- ranger::ranger(Grade ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  pred_Grade <- apply(pred_ranger$prediction,1,which.max)
  temp <- cbind.data.frame(Grade=da_test2$Grade, pred_Grade=pred_Grade, 
                           pred_ranger$prediction, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  colnames(temp)[3:5] <- paste0('Grade',1:3,'_predict_prob')
  res <- rbind.data.frame(res, temp)
  
  res_vus <- VUS(x=temp[temp$Grade==1,5], y=temp[temp$Grade==2,5],
                 z=temp[temp$Grade==3,5], method='NonPar',
                 lam.minus=grade_proportion[1], lam0=grade_proportion[2],
                 lam.plus=grade_proportion[3])
  res_stat <- c(res_stat, res_vus$estimate)
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(VUS=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_Grade.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_Grade.csv'))




###########
## ComBat##
###########

library(sva)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1,check.names = F)
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
saveRDS(da_embedding, file=paste0(input,'embedding_combat.rds'))



########################################
## classification analysis for ComBat ##
########################################

## ER ##
da <- readRDS(paste0(input,'embedding_combat.rds'))
ds <- names(table(meta$DataSet))
test_geo <- paste0('GSE',c(2034, 2990, 3494, 7390))
res_roc <- NULL
res <- res_stat <- c()
for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$ER),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$ER),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  set.seed(123)
  res.boruta <- Boruta(da_train, meta_train$ER)
  
  da_train2 <- cbind.data.frame(ER=factor(meta_train$ER), 
                                da_train[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  da_test2 <- cbind.data.frame(ER=factor(meta_test$ER), 
                               da_test[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  set.seed(123)
  train_ranger <- ranger::ranger(ER ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  er_prob <- pred_ranger$prediction[,1]
  pred_ER <- rep(0,nrow(da_test2))
  pred_ER[er_prob < 0.5] <- 1
  temp <- cbind.data.frame(ER=da_test2$ER, pred_ER=pred_ER, 
                           pred_p=er_prob, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  res <- rbind.data.frame(res, temp)
  
  res_roc[[i]] <- evalmod(scores = 1-temp$pred_p, labels = temp$ER)
  df <- cbind.data.frame(Reference=factor(temp$ER), Prediction=factor(temp$pred_ER))
  res_stat <- c(res_stat, auc(res_roc[[i]])[1,4])
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(AUC=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_ER_ComBat.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_ER_ComBat.csv'))


## Grade ##
test_geo <- paste0('GSE',c(11121, 1456, 2990, 3494, 7390, 9195))
res <- res_stat <- c()
for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Grade),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Grade),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  set.seed(123)
  res.boruta <- Boruta(da_train, meta_train$Grade)
  
  da_train2 <- cbind.data.frame(Grade=factor(meta_train$Grade), 
                                da_train[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  da_test2 <- cbind.data.frame(Grade=factor(meta_test$Grade), 
                               da_test[,names(res.boruta[[1]][res.boruta[[1]]=='Confirmed'])])
  set.seed(123)
  train_ranger <- ranger::ranger(Grade ~ ., data = da_train2, probability = TRUE)
  pred_ranger <- predict(train_ranger, da_test2)
  pred_Grade <- apply(pred_ranger$prediction,1,which.max)
  temp <- cbind.data.frame(Grade=da_test2$Grade, pred_Grade=pred_Grade, 
                           pred_ranger$prediction, dataset=test_geo[i])
  rownames(temp) <- rownames(da_test2)
  colnames(temp)[3:5] <- paste0('Grade',1:3,'_predict_prob')
  res <- rbind.data.frame(res, temp)
  
  res_vus <- VUS(x=temp[temp$Grade==1,5], y=temp[temp$Grade==2,5],
                 z=temp[temp$Grade==3,5], method='NonPar',
                 lam.minus=grade_proportion[1], lam0=grade_proportion[2],
                 lam.plus=grade_proportion[3])
  res_stat <- c(res_stat, res_vus$estimate)
}
res_stat <- formatC(res_stat, digit=3)
res_stat <- cbind.data.frame(VUS=res_stat, DataSet=test_geo)
write.csv(res_stat, file=paste0(input,'stat_summary_Grade_ComBat.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_Grade_ComBat.csv'))




