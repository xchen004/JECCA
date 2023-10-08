# Joint Embedding Canonical Correlation Analysis
JECCA is a method to extract shared information across multiple data sets. The main function is JECCA.R, which can be used by source('JECCA.R'). JECCA_pathway is a function to perform pathway analysis. The script folder includes all codes of the applications of JECCA to three case studies (Breast Cancer, IBD, and RA). The result folder includes all data sets and analysis output for the four case studies.
## Description
### ```JECCA```
JECCA(da,metadata,ncomp=NULL)
#### Arguments

| Arguments   | Description                                                                                                       |
|-------------|-------------------------------------------------------------------------------------------------------------------|
| da          | gene expression matrix of combined multiple datasets. Each column represents a sample, each row represents a gene.|
| metadata    | meta file with first column of sample id matched with column names of gene expression matrix and the second column indicating the dataset information for each sample.|
| ncomp       | total number of JECCA embeddings need to be returned. The number of embeddings should not exceed the minimum number of samples of the multiple datasets.|

#### Value
| Value       | Description                                                                                                       |
|-------------|-------------------------------------------------------------------------------------------------------------------|
| res         | JECCA derived embeddings|


### ```JECCA_pathway```
JECCA_pathway(da,metadata,embedding,pathway_database,output)
#### Arguments

| Arguments   | Description                                                                                                       |
|-------------|-------------------------------------------------------------------------------------------------------------------|
| da          | gene expression matrix of combined multiple datasets. Each column represents a sample, each row represents a gene.|
| metadata    | meta file with first column of sample id matched with column names of gene expression matrix and the second column indicating the dataset information for each sample.|
| embedding   | JECCA output embeddings.|
| pathway_database   | The pathway database need to be selected for the GSEA pathways analysis. The following pathway databases are provided: ReactomPA, KEGG, GOALL, GOBP, GOMF, GOALL. User can select any combinations of the above databases to perform GSEA pathway analysis. fgsea function is used to perform reactomPA GSEA analysis. gseKEGG and gseGO function is used to perform KEGG and GO GSEA analysis.|
| output      | The output directory. The function will generate following results in the output directory: subtype expression based on the da matrix and embedding, GSEA pathway results.|
