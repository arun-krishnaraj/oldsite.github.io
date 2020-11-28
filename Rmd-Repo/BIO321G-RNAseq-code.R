source("bio321g_rnaseq_utils.R")

library(DESeq2)
## DESeq object creation and testing
dds = DESeqDataSetFromMatrix(countData = rnaCounts, colData = sampleAnnotation, 
                             design = ~ time + genotype + time:genotype)

dds = DESeq(dds, parallel = FALSE, test = "LRT", reduced = ~ time + genotype)
ddsresults <- results(dds, alpha = 0.1) 

## FDR calculation
which(!is.na(ddsresults$padj)) -> testedGenes
sum(ddsresults[testedGenes, "padj"] <= .1)

## Count normalization and transformation 
lgNorm = log2(counts(dds, normalized = TRUE)+1)

## PCA on log normalized counts, all genes 
pca = prcomp(t(lgNorm))

pca.data <- data.frame(pca$x[,1:2])
pca.data$group  = sampleAnnotation[rownames(pca.data), "group"]
pca.data$sample = rownames(pca.data)

pca$sdev[[1]]^2/sum(pca$sdev^2) # Fraction of variance calculations
pca$sdev[[2]]^2/sum(pca$sdev^2)

library(ggplot2); library(tidyverse)
pca1 <- pca.data %>% ggplot(aes(x = PC1, y = PC2, color = group, label = sample)) + 
  geom_point() + scale_color_manual(values=groupColors) + theme_minimal() + theme(panel.background = element_rect(fill = NULL, color = "black")) +
  ggtitle("PCA Plot for all Genes") + labs(color = "Group") + xlab("PC1 (55.9% explained var.)") + ylab("PC2 (15.5% explained var.)")
pca1

## Subset to include assigned genes
geneset <- read.table("gene_sets.tsv.gz", 
                      sep="\t", row.names=NULL, header=TRUE, quote="", comment.char = "")
gene_vec <- geneset[geneset$gene_ontology_primary_id == "GO:0006865",]$gene 

## Assigned geneset information, tsv creation
gene_df <- geneNamesAndDescriptions[geneNamesAndDescriptions$gene %in% gene_vec, 1:3] %>% data.frame() 
rownames(gene_df) <- NULL

write.table(gene_df, "Final-Project-Gene-desc.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

## Filtering counts for assigned geneset
lgGo <- lgNorm[rownames(lgNorm) %in% gene_vec,] %>% data.frame(check.names = FALSE) 

## PCA on log normalized counts, assigned genes
pca_lggo <- prcomp(t(lgGo))

pca_lggo_dat  <- data.frame(pca_lggo$x[,1:2])
pca_lggo_dat$group <- sampleAnnotation[rownames(pca_lggo_dat), "group"]
pca_lggo_dat$sample <- rownames(pca_lggo_dat)

pca_lggo$sdev[[1]]^2 / sum(pca_lggo$sdev^2) # Fraction of variance calculations
pca_lggo$sdev[[2]]^2 / sum(pca_lggo$sdev^2)

pca2 <- pca_lggo_dat %>% ggplot(aes(x = PC1, y = PC2, color = group, label = sample)) + 
  geom_point() + scale_color_manual(values=groupColors) + ggtitle("PCA Plot for Gene Set (GO:0006865)") + 
  labs(color = "Group") + theme_minimal() + theme(panel.background = element_rect(fill = "white", color = "black")) +
  xlab("PC1 (54.9% explained var.)") + ylab("PC2 (22.1% explained var)")
pca2

## Clustered heatmap, assigned genes
library(pheatmap)
heatData <- lgGo -rowMeans(lgGo)
heatData[heatData > 2] = 2
heatData[heatData < -2] = -2

pheat1 <- pheatmap(heatData, color = heatPalette, clustering_method = "average", 
         labels_row = geneNamesAndDescriptions[rownames(heatData), "symbol"], 
         cutree_cols = 2, fontsize_row = 8, main = "Clustered Heatmap") 
pheat1

## Expression stripchart, assigned genes 
data.frame(results(dds)[5]) -> pmat

top_9_exp <- subset(pmat, rownames(pmat) %in% gene_vec) %>% arrange(pvalue) %>% head(9) %>% rownames 

GoExp <- lgGo[rownames(lgGo) %in% top_9_exp,]
stripgg <- stripchart321g(data = GoExp, sampleAnnotation = sampleAnnotation)
stripgg <- stripgg + ggtitle("Gene Expression by Genotype and Time") + xlab("Genotype") + ylab("Expression") + labs(color = "Group", shape = "Time")
  
  
  
  