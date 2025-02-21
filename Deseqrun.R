#scripts to run differential gene expression RNA seq analysis using Deseq package 
#loading the required libraries
library(utils)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages("gsub")
library(DESeq2)
library(magrittr)
library(gsub)

#importing and unzipping the data
counts<- read.table("C:/Users/olabo/PRACTICE_BINF/Graduate_Project_RNASeq/Bioinformatics-_Final-_Project/GSE260666_raw_counts_GRCh38.p13_NCBI.tsv.gz", 
                        header= TRUE,
                        sep ="\t",
                        row.names=1)
metadata<- read.csv("metadata.csv")
metadata<- metadata[-7, ] #GSM8122041 is missing from the raw counts data i downloaded and if the column to rowname of the count to metadata does not match  deseq cannot run 


head(counts)
head(metadata)
#Wrangle data fOr desq 

counts$GeneID <- rownames(counts)  #ConvertS row names to a column 
counts <- counts[, c(ncol(counts), 1:(ncol(counts) - 1))]  # MOVES GENE ID TO THE FIRAST COLUN 
colnames(counts)[1] <- "GeneID"  # Rename columns (first column as "GeneID")

#Deseq needs the counts to have gene iDs as row names 
rownames(counts)=counts$GeneID  # assigns geneid as the row names 
counts2= counts[,-c(1,2)] # removes column one and 2 from the data as deseq =2 does not require the information 


#sample id as rowname for also metadata
rownames(metadata)<- metadata$SampleID # sets sampleid as row name for metadata

metadata<- metadata[, -1, drop = FALSE]# Drop included so it maintains a data frame format 


# removing space from meta data as it affects deseq and to aviod warnings 
metadata$Condition[metadata$Condition== "Control_rep1 [21120D1]"]= "Control_rep1"
metadata$Condition <- gsub(" \\[.*\\]", "", metadata$Condition)  #  regular expression (regex), which helps find and replace specific patterns in text.
#" " (Space), [ is a special character in regex, so we escape it with \\. .* means "match everything after the [ until the next part of the pattern".


# Factor conditions in the metadata for deseq analysis 
metadata$Condition <- gsub("_rep[0-9]+", "", metadata$Condition) #_rep[0-9]+ → Finds "_rep" followed by any number (removes it).
metadata$Condition <- factor(metadata$Condition)
levels(metadata$Condition)


#Deseq Run 

dds<-DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)

