#scripts to run differential gene expression RNA seq analysis using Deseq package 
#loading the required libraries
library(utils)
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #install.packages("BiocManager")}
install.packages("tidyverse")
BiocManager::install("EnhancedVolcano")
install.packages("EnhancedVolcano")
install.packages("msigbr")

BiocManager::install("enrichplot")
BiocManager::install("GSAR")
library(fgsea)
library(apeglm)
library(DESeq2)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(GSAR)
library(enrichplot)


 #importing and unzipping the data
counts <- read.delim("C:/Users/olabo/PRACTICE_BINF/Graduate_Project_RNASeq/Bioinformatics-_Final-_Project/GSE260666_raw_counts_GRCh38.p13_NCBI.tsv.gz", 
                     header = TRUE,
                     sep= "\t",
                     row.names = 1)

metadata<- read.csv("metadata.csv")
metadata<- metadata[-7, ] #GSM8122041 is missing from the raw counts data i downloaded and if the column to rowname of the count to metadata does not match  deseq cannot run 

head(counts)
head(metadata)
counts2<- counts

#Wrangle data fOr desq 
counts2$GeneID <- rownames(counts2)  # Move GeneID to a column if not already
counts2 <- counts2[ c(ncol(counts2), 1:(ncol(counts2) - 1))]  # Move GeneID column to first
colnames(counts2)[1] <- "GeneID"  # Rename the first column as "GeneID"

# Now, set GeneID as row names and remove it as a column
#rownames(counts) <- counts$GeneID  # Set GeneID as the row names
#counts <- counts[, -1]  # removes column one and 2 from the data as deseq =2 does not require the information 


#Getting the gene names of the gene _ids used in this study from ensemble 

#counts2 will be used to retrieve the gene names from ensemble
counts3<- counts2%>% 
  select(GeneID)


rownames(counts3)<-NULL #removing the row names

 

# writing out the file 
write.csv(genes_names, "genenames.csv")


#sample id as row name for also metadata
rownames(metadata)<- metadata$SampleID # sets sample-id as row name for metadata

metadata<- metadata[, -1, drop = FALSE]# Drop included so it maintains a data frame format 


# removing space from metadata as it affects deseq and to avoid warnings 
metadata$Condition[metadata$Condition== "Control_rep1 [21120D1]"]= "Control_rep1"
metadata$Condition <- gsub(" \\[.*\\]", "", metadata$Condition)  #  regular expression (regex), which helps find and replace specific patterns in text.
#" " (Space), [ is a special character in regex, so we escape it with \\. .* means "match everything after the [ until the next part of the pattern".


# Factor conditions in the metadata for deseq analysis 
metadata$Condition <- gsub("_rep[0-9]+", "", metadata$Condition) #_rep[0-9]+ → Finds "_rep" followed by any number (removes it).
metadata$Condition <- factor(metadata$Condition)
levels(metadata$Condition)



#Deseq Run 



dds<-DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)

dds<- dds[rowSums(counts(dds)) > 10]  #ignores gennes with counts lower than 10 
dds<-DESeq(dds)  #Running deseq

res<- results(dds)

res2<- res #making a duplicate of res

resultsNames(dds)



# comparing NAfL VS Control
result_NAFL= results(dds, contrast = c("Condition", "NAFL", "Control"), alpha = 0.05)


summary(result_NAHL)

result_NAFL["6566",]

# Compare NASH vs. Control
results_NASH = results(dds, contrast = c("Condition", "NASH", "Control"), alpha = 0.05)
summary(results_NASH)

# Compare NASH vs. NAFL
results_NASH_NAFL = results(dds, contrast = c("Condition", "NASH", "NAFL"), alpha = 0.05)
summary(results_NASH_NAFL)

# GETTING GENE NAMES USING ANOTATIONdb

columns(org.Hs.eg.db)
anno<- AnnotationDbi::select(org.Hs.eg.db, rownames(res),columns = c("ENSEMBL","ENTREZID","SYMBOL","GENENAME"),
                             keytype= "ENTREZID")

anno_unique <- anno[!duplicated(anno$ENTREZID), ]

#res_annotated<- merge(as.data.frame(res), anno_unique, by.x= "row.names", by.y= "ENTREZID")

res2_df <- as.data.frame(res)
res2_df$ENTREZID<- rownames(res2_df)
res_annotated<- merge(res2_df, anno_unique,by= "ENTREZID")



#-- Getting the overexpressed gene 
sig_gene<- res_annotated %>% 
  dplyr::filter(padj<0.05& log2FoldChange>2) %>% 
  arrange(padj) %>% 
  na.omit() %>%
  select("GENENAME","ENTREZID", "SYMBOL") 

#--- Get under-expressed genes 
unde_gene<-res_annotated %>% 
  dplyr::filter(padj<0.05 & log2FoldChange <2) %>% 
  arrange(padj) %>% 
  na.omit() %>% 
  select("GENENAME","ENTREZID")
  







#queries for scfa genes

#subset(res_annotated, grepl("^HDAC3$", SYMBOL))
scfa_genes<- c("FFAR2", "FFAR3", "SLC16A1", "HDAC3","SLC16A3", "SLC5A8", "ABCG2",
               "IL22", "SLC5A8", "HCAR2", "HDAC2", "PPARG", "NFKB1", "IL10", "SLC16A4", "TJP1","CLDN1")

scfa_genes_deseq_res<-subset(res_annotated, SYMBOL %in% c("FFAR2", "FFAR3", "SLC16A1", "HDAC3","SLC16A3", "SLC5A8", "ABCG2",
                                                          "IL22", "SLC5A8", "HCAR2", "HDAC2", "PPARG", "NFKB1", "IL10", "TNF", "SLC16A4", "TJP1","CLDN1"))
#ffar3 Count 
                                                          
plotCounts(dds, gene="2865", intgroup="Condition")    

# ffar2 count 

plotCounts(dds, gene="2867", intgroup = "Condition")

# enhanced volcano plot for 
library(EnhancedVolcano)

EnhancedVolcano(res_annotated, lab = res_annotated$SYMBOL, 
                x = "log2FoldChange", y = "pvalue", 
                selectLab = scfa_genes,  # Highlight SCFA genes
                pCutoff = 0.05, FCcutoff = 1)

ggsave("Volcano_plot_SCFA_genes.png")



# Extract expression values for HDAC3 (gene ID 8841)
hdac3_expr <- counts["8841", ]

# Create a data frame for plotting
hdac3_df <- data.frame(
  Expression = as.numeric(hdac3_expr),
  Condition = metadata$Condition
)

# Load ggplot2
library(ggplot2)

# Plot boxplot
ggplot(hdac3_df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "HDAC3 Expression Across Conditions",
       y = "Normalized Expression",
       x = "") +
  theme_minimal()



#Pathway analysis

# get Hallmark gene sets

gene_sets<-msigdbr(specie="human", category= "C2")

gene_sets<- gene_sets %>% 
  dplyr::select(gs_name, gene_symbol)

#-- Getting the over expressed gene 
sig_gene<- res_annotated %>% 
  dplyr::filter(padj<0.05& log2FoldChange>2) %>% 
  arrange(padj) %>% 
  dplyr::pull(SYMBOL) %>% 
  na.omit()

#--- Get under-expressed genes 
res_annotated %>% 
  dplyr::filter(padj<0.05 & log2FoldChange <2) %>% 
  write.csv(file = "under_expressed_geness.csv")


# over-representation analysis using enricher 
egmt<- enricher(gene= sig_gene,
                TERM2GENE = gene_sets)


edf<- as.data.frame(egmt)
view(edf)

#plot result with cluster 
dotplot(egmt, title = "Over-Represented Pathways (ORA)") +
  theme(axis.text.y = element_text(size = 8))
barplot(egmt, title = "Over-Represented Pathways (ORA)") +
  theme(axis.text.y = element_text(size = 8))




#GSEA

gsea_res <- res_annotated %>% 
  arrange(padj) %>% 
  mutate(gsea_metric = sign(log2FoldChange)* -log10(padj)) %>%
  arrange(desc(gsea_metric)) %>% 
  select(SYMBOL, gsea_metric) %>% 
  na.omit()



#turns it into a named vector, with: names = SYMBOLs values = gsea_metric
rank<- gsea_res %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  deframe()

head(rank)

# run GSEA
gseares<- GSEA(geneList = rank,
               TERM2GENE = gene_sets )

gseares_df<- as.data.frame(gseares)

class(gsea_res)

View(gseares_df)

#summary plots
dotplot(gseares)

# gsea plot for 

gseaplot(gseares, geneSetID = "REACTOME_TRANSLATION",
         title= "REACTOME_TRANSLATIONPathway ")

# gsea plot for top overrepresented genes  

top_pathway<- gseares_df%>% 
  top_n(n=4, wt=NES) %>% 
  pull(ID)
head(top_pathway)

# make gseapot for each and return a list 
top_pathway_plots<- lapply(top_pathway, function(pathway){
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})


#--- arrange with labels as multi- panel plot 
top_pathway_plot<- ggarrange(plotlist= top_pathway_plots, ncol=2, nrow=2, labels= "AUTO")
print(top_pathway_plot)
top_pathway_plot

library(ggpubr)
library(ggarrange)

gseares<- GSEA(geneList = rank,
               TERM2GENE = scfa_pathways )
fgseaRes <- fgsea(pathways = scfa_pathways, stats = rank)

#Gsar analysis



library(GSAR)
scfa_genes <- list(
  SCFA_Pathway = c("GPR43", "GPR41", "HDAC9", "ACADS", "ACADL", "SLC5A8", "ACOT1")
)

# Rank genes by log2FC
ranked_genes <- sort(res_annotated$log2FoldChange, decreasing = TRUE)
# Run GSAR with custom SCFA gene set
names(ranked_genes)<- (res_annotated$SYMBOL)
gsar_res <- TestGeneSets(ranked_genes, geneset = scfa_genes)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
gsar_res <- TestGeneSets(ranked_genes, scfa_genes)


library(clusterProfiler)
library(msigdbr)

# Rank genes by log2FC
gene_ranks <- res_annotated$log2FoldChange
names(gene_ranks) <- res_annotated$GENENAME

gene_ranks <- sort(gene_ranks, decreasing = TRUE)
gene_ranks <- tapply(gene_ranks, names(gene_ranks), mean)
gene_ranks <- gene_ranks[!is.na(names(gene_ranks)) & !is.na(gene_ranks)]
gene_ranks <- as.numeric(gene_ranks)
names(gene_ranks) <- names(gene_ranks)
# Get SCFA-related gene sets (KEGG/Reactome)
scfa_setss <- msigdbr(species = "human", category = "C2") %>% 
  filter(grepl("butanoate|PPAR|fatty_acid|HDAC", gs_name, ignore.case = TRUE))

# Run GSEA
gsea_res <- GSEA(rank, TERM2GENE = scfa_setss[, c("gs_name", "gene_symbol")], pvalueCutoff = 0.1)
gsea_res_df<- as.data.frame(gsea_res) 

[1] "SENESE_HDAC1_TARGETS_UP"                     
[2] "PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP"
[3] "PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_DN"
[4] "LUI_TARGETS_OF_PAX8_PPARG_FUSION"
dotplot(gsea_res, showCategory = 10) + ggtitle("GSEA: SCFA-Related Pathways")
ggsave("GSEA_SCFA_pathways.png")

gsea_resdf<-as.data.frame(gsea_res)
gseaplot(gsea_res, geneSetID = "SENESE_HDAC1_TARGETS_UP",
         title= "SENESE_HDAC1_TARGETS_UP Pathway ")

# gsea plot for top overrepresented genes  

top_pathwayy<- gsea_resdf%>% 
  top_n(n=4, wt=NES) %>% 
  pull(ID)
head(top_pathwayy)



scfa_sets <- list(
  # Receptors (G-protein-coupled receptors)
  Receptors = c("FFAR2", "FFAR3", "HCAR2", "GPR109A", "GPR41", "GPR43"),
  
  # Transporters (monocarboxylate/SLC family)
  Transporters = c("SLC16A1", "SLC16A3", "SLC16A4", "SLC5A8", "ABCG2"),
  
  # Signaling (downstream mediators)
  Signaling = c("HDAC3", "HDAC1", "PPARG", "NFKB1", "IL10", "IL22", "TNF")
)



#Log fold change shrinkage for visualization and ranking
#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
# "Condition_NAFL_vs_Control",
resLFC<- lfcShrink(dds, coef = "Condition_NASH_vs_Control", type = "apeglm")
reslfc_naflVScontrol<- lfcShrink(dds, coef = "Condition_NAFL_vs_Control", type = "apeglm"  )

# MA plot 
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, ylim=c(-2,2))

plotMA(reslfc_naflVScontrol,ylim=c(-2,2))



# Visualization 
## 1 dispersion plot, as a genes read count increases  dispersion decreases.

plotDispEsts(dds)


#PCA 
# explains variance transformation 


# Variance Stabilizing transformation
gene_symbols<- res_annotated$SYMBOL
# Make a named vector mapping ENTREZID to SYMBOL

names(gene_symbols) <- res_annotated$ENTREZID
# Rename the rownames of vsd using the SYMBOL names
rownames(vsd) <- gene_symbols[rownames(vsd)]
                              
sample_annots<-as.data.frame(metadata$Condition) 
rownames(sample_annots) <- colnames(vsd) 

vsd <- vst(dds, blind= FALSE)

sub_count<- assay(vsd)[rownames(vsd)%in% scfa_genes,]
head(rownames(vsd))

# heat map for scfa gene expression across samples 
pheatmap(sub_count,
         annotation_col = sample_annots,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "SCFA Gene Expression Across Samples")


#transformed values are used to generate a PCA  plot 

plotPCA(vsd, intgroup= "Condition")

#using rdlogs 
rld<- rlog(dds, blind = FALSE)
pca_data<- plotPCA(rld, intgroup= c("Condition")) #returnData= TRUE
plot(pca_data$x[,3], pca_data$x[,4], col=Condition, pch=16)

plotPCA(rld, intgroup= c("Condition"))
#heatmaps.

### Heat maps of sample to sample using a distance matrix with clustering ###

#distance matrix 
sampledist<- dist(t(assay(vsd)))
sampledist_matrix<- as.matrix(sampledist)
colnames(sampledist_matrix)
#

# setting a color scheme
color<- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)


# generating the heatmap of sample to sample 
pheatmap(sampledist_matrix, clustering_distance_rows = sampledist,
         clustering_distance_cols = sampledist, col = color)
# notes lighter colors = less similar 

# volcano plot  for degs
EnhancedVolcano(res_annotated, 
                lab = gene_symbols,
                x= "log2FoldChange",
                y= "padj",
                pCutoff = 0.05,
                FCcutoff =1,
                pointSize = 2,
                labSize = 4,
                title = "Differential Expresion in NAFLD ",
                subtitle = 'padj < 0.05 & |log2FC| > 1',
                caption = 'DESeq2 results: GSE260666',
                colAlpha = 1,
                selectLab = c("HDAC3", "FFAR2", "FFAR3", "MCT1", "MCT4", "SMCT1", "ABCG2"))
