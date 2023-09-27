# Joint Embedding Canonical Correlation Analysis
JECCA is a method to extract shared information across multiple data sets. The main function is JECCA.R, which can be used by source('JECCA.R'). JECCA_pathway is a function to perform pathway analysis. The script folder includes all codes of the applications of JECCA to three case studies (Breast Cancer, IBD, and RA). The data folder includes all data sets for the three case studies.
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
