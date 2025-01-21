#scripts to run differential gene expression RNA seq analysis using Deseq package 
#loading the required libraries
library(utils)
install.packages('DESEQ2')

#importing and unzipping the data
tpm_counts<- read.table("C:/Users/olabo/PRACTICE_BINF/Graduate_Project_RNASeq/Bioinformatics-_Final-_Project/GSE260666_norm_counts_TPM_GRCh38.p13_NCBI.tsv", 
                        header= TRUE,
                        sep ="\t",
                        row.names=1)
head(tpm_counts)
