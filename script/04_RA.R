module load hdf5
module load gcc/10.2.0
module load zlib
export PATH=$PATH:/projects/abv/DSS/xchen/app/R/R-4.3.0/bin/
#export PATH=$PATH:/sc/wo/home/chenxx25/app/R/R-4.1.2/bin/
R


library(plyr)
library(dplyr)
library(reshape2)
library(GEOquery)
library(Seurat)


####################################
## data loading and preprocessing ##
####################################

input <- '/.../Data/RA/'
GEOs <- c("GSE129705","GSE138746_ada_pbmc","GSE33377","GSE42296",
          "GSE54629", "GSE12051")

phenotData <- c()
#Load gset and annot files to the environment
for(i in GEOs){
  load(paste0(input,i,".gset"))
  assign(paste0("data",i), da)
  load(paste0(input,i,".annot")) 
  anno <- cbind.data.frame(sid=anno$sid, trt=anno$trt, 
                           Group=anno$Group, GEO=i)
  phenotData <- rbind(phenotData, anno)
  print(sum(is.na(da)))}

#Merge datasets
genes_in_common <- c()
for(i in GEOs){
  if(length(genes_in_common)==0){ genes_in_common <- rownames(get(paste0("data",i)))
  } else { genes_in_common <- genes_in_common[genes_in_common %in% rownames(get(paste0("data",i)))] } }
data <- matrix()
for(i in GEOs){
  if(dim(data)[1] == 1){ data <- get(paste0("data",i))[genes_in_common,]
  } else { data <- cbind(data,get(paste0("data",i))[genes_in_common,]) }}

rownames(phenotData) <- phenotData$sid
phenotData <- phenotData[,c(4,3)]
data <- data[,rownames(phenotData)]
levels(phenotData$Group) <- c('NR', 'R', 'R')

da_s <- anno_all <- da_all <- da_all2 <- NULL
for(i in 1:length(GEOs))
{
  temp_anno <- phenotData[phenotData$GEO %in% GEOs[i],]
  temp_da <- data[,rownames(temp_anno)]
  da_s[[i]] <- t(scale(t(temp_da)))
  anno_all <- rbind(anno_all, temp_anno)
  da_all <- rbind(da_all, t(temp_da))
  da_all2 <- rbind(da_all2, t(da_s[[i]]))
}
anno_all <- cbind.data.frame(sid=rownames(anno_all),anno_all)
colnames(anno_all)[2] <- 'DataSet'
output <- '/.../Results/RA/'
dir.create(output, recursive=T)
write.csv(anno_all, file=paste0(output,'anno_hg.csv'),row.names=F)
write.csv(t(da_all), file=paste0(output,'exprs_hg.csv'))
write.csv(t(da_all2), file=paste0(output,'exprs_scale_hg.csv'))



##############
## multiCCA ##
##############
source('/.../Code/JECCA.R')
input <- '/.../Results/RA/'
da <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1,check.names = F)
meta <- read.csv(paste0(input,'anno_hg.csv'))
out <- JECCA(as.matrix(da),meta,ncomp=20)
saveRDS(out, file=paste0(input,'embedding.rds'))


## t-test between R/NR and CCA1-10 ##

da_embedding_all <- readRDS(paste0(input,'embedding.rds'))
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
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


#############################
## classification analysis ##
#############################

library(precrec)
ds <- names(table(meta$DataSet))
test_geo <- ds
da <- readRDS(paste0(input,'embedding.rds'))
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:5])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:5])
  
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
write.csv(res_stat, file=paste0(input,'stat_summary.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred.csv'))


###########
## ComBat##
###########

library(sva)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1,check.names = F)
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

library(precrec)
ds <- names(table(meta$DataSet))
test_geo <- ds
da <- readRDS(paste0(input,'embedding_combat.rds'))
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:5])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:5])
  
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



############
## QN ##
############

library(preprocessCore)
meta <- read.csv(paste0(input,'anno_hg.csv'),row.names = 1)
da_all <- read.csv(paste0(input,'exprs_hg.csv'),row.names=1,check.names = F)
da_all <- as.matrix(da_all)
GEOs <- GEOs[c(6,1:5)]
meta2 <- meta[meta$DataSet %in% GEOs[1],]
da_QN <- normalize.quantiles(da_all[,rownames(meta2)])
rownames(da_QN) <- rownames(da_all)
colnames(da_QN) <- rownames(meta2)
qn.ref <- normalize.quantiles.determine.target(da_all[,rownames(meta2)],
                                               target.length = nrow(da_all))
for(i in 2:length(GEOs))
{
  temp_meta <- meta[meta$DataSet %in% GEOs[i],]
  temp_da <- da_all[,rownames(temp_meta)]
  temp_QN <- normalize.quantiles.use.target(temp_da, qn.ref, copy = F)
  da_QN <- cbind.data.frame(da_QN,temp_QN)
}
da_QN <- da_QN[,colnames(da_all)]

da_seurat <- CreateSeuratObject(counts=da_QN)
hvgs <- rownames(da_seurat)
da_seurat <- ScaleData(da_seurat)
da_seurat <-  RunPCA(da_seurat, features = hvgs, npcs = 30)
da_seurat@meta.data <- cbind.data.frame(da_seurat@meta.data, 
                                        anno_all[match(colnames(da_seurat), rownames(anno_all)),])

da_seurat <- da_seurat %>% RunUMAP(reduction = "pca",dims=1:30) %>% 
  FindNeighbors(reduction = "pca",graph.name="rna.snn.h",dims=1:30)
da_embedding <- Embeddings(da_seurat, reduction = 'pca')
saveRDS(da_embedding, file=paste0(input,'embedding_qn.rds'))


########################################
## classification analysis for QN ##
########################################

library(precrec)
ds <- names(table(meta$DataSet))
test_geo <- ds
da <- readRDS(paste0(input,'embedding_qn.rds'))
res_roc <- NULL
res <- res_stat <- c()

for(i in 1:length(test_geo))
{
  meta_train <- meta[meta$DataSet != test_geo[i] & !is.na(meta$Group),]
  meta_test <- meta[meta$DataSet == test_geo[i] & !is.na(meta$Group),]
  da_train <- da[rownames(meta_train),]
  da_test <- da[rownames(meta_test),]
  da_train2 <- cbind.data.frame(Group=factor(meta_train$Group), da_train[,1:5])
  da_test2 <- cbind.data.frame(Group=factor(meta_test$Group), da_test[,1:5])
  
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
write.csv(res_stat, file=paste0(input,'stat_summary_qn.csv'),row.names=F)
write.csv(res, file=paste0(input,'res_pred_qn.csv'))




#####################################
## xCell cell composition analysis ##
#####################################

library(xCell)
input <- '/.../Data/RA/'
GEOs <- c("GSE129705","GSE138746_ada_pbmc","GSE33377","GSE42296",
          "GSE54629", "GSE12051")
#Load gset and annot files to the environment
for(i in GEOs){
  load(paste0(input,i,".gset"))
  
  data <- da[!duplicated(rownames(da)),]
  if(i == 'GSE129705')
  {res <- xCellAnalysis(data)}
  else
  {res <- cbind.data.frame(res,xCellAnalysis(data))}
}

res <- res[,rownames(meta)]
output <- '/.../Results/RA/'
write.csv(res, file=paste0(output,'xcell.csv'))
