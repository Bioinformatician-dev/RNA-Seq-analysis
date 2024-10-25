# Step 1: Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("ggplot2", "pheatmap"))

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  # For human data

# Step 2: Simulate Example Data
set.seed(42)  # For reproducibility
n_genes <- 1000
n_samples <- 6

# Create a count matrix
counts <- matrix(rpois(n_genes * n_samples, lambda = 10), nrow = n_genes)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Sample", 1:n_samples)

# Create sample metadata
colData <- data.frame(
  condition = factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")),
  row.names = colnames(counts)
)

# Step 3: Create DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# Step 4: Preprocessing and Normalization
dds <- DESeq(dds)

# Step 5: Differential Expression Analysis
res <- results(dds)
summary(res)

# Step 6: Visualize Results

## MA Plot
ggplot(as.data.frame(res), aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  xlab("log2(Mean Expression)") +
  ylab("log2(Fold Change)") +
  theme_minimal() +
  ggtitle("MA Plot")

## Volcano Plot
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  xlab("log2(Fold Change)") +
  ylab("-log10(Adjusted p-value)") +
  theme_minimal() +
  ggtitle("Volcano Plot")

# Step 7: Heatmap of Significant Genes
sig_genes <- res[which(res$padj < 0.05), ]
mat <- assay(dds)[rownames(sig_genes), ]

# Create heatmap
pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Heatmap of Significant Genes")

# Step 8: Functional Enrichment Analysis
gene_list <- rownames(sig_genes)

# Perform GO enrichment analysis
go_results <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", 
                        ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# View GO results
print(head(go_results))

# Plot GO results
barplot(go_results, showCategory = 10, title = "GO Enrichment Analysis")

# Step 9: Save Your Results
write.csv(as.data.frame(res), file = "differential_expression_results.csv")
write.csv(as.data.frame(sig_genes), file = "significant_genes.csv")

