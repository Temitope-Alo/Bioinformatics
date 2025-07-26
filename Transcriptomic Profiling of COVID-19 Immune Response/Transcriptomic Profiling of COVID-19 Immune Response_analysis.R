# GSE152418
setwd("C:/Users/temit/Desktop/My_GitHub/ID_bioinformatics/1. GSE152418")

# Load required packages
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(EnhancedVolcano)

# Step 2: Load and Prepare the Data
counts <- read.delim(gzfile("GSE152418_p20047_Study1_RawCounts.txt.gz"), 
                     row.names = 1, check.names = FALSE)
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Step 3: Create sample metadata
sample_ids <- colnames(counts)
group <- c(rep("COVID-19", 17), rep("Convalescent", 17))
pheno <- data.frame(sample_id = sample_ids, group = factor(group))
rownames(pheno) <- pheno$sample_id

# Step 4: Match sample order
all(colnames(counts) %in% rownames(pheno))  # Should be TRUE
counts <- counts[, rownames(pheno)]

# Step 5: Create DGEList and filter
dge <- DGEList(counts = counts, group = pheno$group)
keep <- filterByExpr(dge) # Filter lowly expressed genes
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Library Size Barplot: Check sequencing depth before normalization.
barplot(dge$samples$lib.size / 1e6, names.arg = colnames(dge), 
        las = 2, cex.names = 0.7, col = "skyblue",
        main = "Library Sizes", ylab = "Million Reads")

# Step 6: Normalization
dge <- calcNormFactors(dge)

# Boxplot of Normalized Counts:Check distribution across samples after normalization.
logCPM <- cpm(dge, log = TRUE)
boxplot(logCPM, las = 2, col = "lightgreen",
        main = "Log2 CPM After Normalization",
        ylab = "Log2 Counts per Million")

# Multidimensional Scaling (MDS) Plot: Assess similarity between samples based on gene expression.
plotMDS(dge, col = as.numeric(dge$samples$group))
legend("topright", legend = levels(dge$samples$group), 
       col = 1:length(levels(dge$samples$group)), pch = 16)

# Step 7: Design matrix and modeling
design <- model.matrix(~group, data = pheno)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit)

# Mean-Variance Plot: Diagnostic plot after model fitting.
plotQLDisp(fit)

# Step 8: Differential Expression Results
topTags(qlf)
results <- topTags(qlf, n = Inf)
summary(decideTests(qlf))
write.csv(results$table, "DEGs_GSE152418.csv")

# Smear Plot (MA Plot): Visualize fold change and abundance.
plotSmear(dge, de.tags = rownames(results$table)[results$table$FDR < 0.05])
abline(h = c(-1, 1), col = "red")

# Volcano Plot: Assess log fold change vs significance.

# Assuming 'results' is your differential expression results object
# and contains: logFC, PValue, FDR, and gene names as rownames

results_df <- as.data.frame(results)

# Basic EnhancedVolcano plot
EnhancedVolcano(
  results_df,
  lab = rownames(results_df),           # gene names as labels
  x = 'logFC',                          # column with log2 fold changes
  y = 'PValue',                         # column with p-values
  pCutoff = 0.05,                       # FDR threshold (for y-axis significance)
  FCcutoff = 1,                         # Fold-change threshold
  title = "Volcano Plot of DEGs",
  subtitle = "FDR < 0.05 and |logFC| > 1",
  xlab = "log2 Fold Change",
  ylab = "-log10(p-value)",
  pointSize = 2.0,
  labSize = 3.5,
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,               # Draw lines connecting labels to points
  widthConnectors = 0.5
)   # Highlighting significant DEGs (FDR < 0.05 and |logFC| > 1)


# Heatmap of Top DEGs
group <- factor(dge$samples$group)   # Ensuring group is a factor
top_genes <- rownames(results)[1:50]
mat <- logCPM[top_genes, ]

annotation_col <- data.frame(Group = group)  # Creating annotation_col with proper rownames
rownames(annotation_col) <- colnames(mat)

pheatmap(mat,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Top 50 DEGs Heatmap")   # Top 50 genes, clustered across samples.



# ----------------------------------------------------------
# Functional Enrichment (GO, KEGG)

# Extract Significant Genes: Let’s filter genes with FDR < 0.05 and |logFC| > 1:
sig_genes <- subset(results_df, FDR < 0.05 & abs(logFC) > 1)
gene_symbols <- rownames(sig_genes)

# Convert Gene Symbols to Entrez IDs
gene_entrez <- bitr(gene_symbols, 
                    fromType = "ENSEMBL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)


# Check the mapped genes:
head(gene_entrez)

# Gene Ontology (GO) Enrichment
go_enrich <- enrichGO(gene         = gene_entrez$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      ont          = "BP",  # BP = Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

head(go_enrich)  # View results

# Visualize GO
barplot(go_enrich, showCategory = 15, title = "GO Biological Process")
dotplot(go_enrich, showCategory = 15, title = "GO Enrichment")

# KEGG Pathway Enrichment
kegg_enrich <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = 'hsa',  # Human
                          pvalueCutoff = 0.05)

head(kegg_enrich)  # View results

# Visualize KEGG
barplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")
dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Dotplot")

# Save Results
write.csv(as.data.frame(go_enrich), "GO_enrichment_results.csv")
write.csv(as.data.frame(kegg_enrich), "KEGG_enrichment_results.csv")


# Enrichment Map (network-like plots)
# Compute similarity between terms
go_enrich_sim <- pairwise_termsim(go_enrich)

emapplot(go_enrich_sim, showCategory = 30, layout = "kk")  # Kamada-Kawai layout
cnetplot(go_enrich, showCategory = 5)



# ----------------------------------------------------------
# GSEA 

# Prepare ranked gene list from DE results
results_df <- as.data.frame(results$table)  # Ensure correct extraction
results_df <- results_df[!is.na(results_df$logFC), ]
gene_list <- results_df$logFC
names(gene_list) <- rownames(results_df)  # Ensembl IDs as names
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert Ensembl IDs to Entrez IDs using BioMart
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_conversion <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = names(gene_list),
  mart = mart
)

# Clean mapping: remove NA and duplicates
gene_conversion <- gene_conversion[!is.na(gene_conversion$entrezgene_id), ]
gene_conversion <- gene_conversion[!duplicated(gene_conversion$ensembl_gene_id), ]

# Filter and rename ranked gene list using Entrez IDs
gene_list_filtered <- gene_list[gene_conversion$ensembl_gene_id]
names(gene_list_filtered) <- gene_conversion$entrezgene_id

# Final clean-up: remove duplicated Entrez IDs (if any)
gene_list_filtered <- gene_list_filtered[!duplicated(names(gene_list_filtered))]

# Final sorting (already sorted, but for safety)
gene_list_filtered <- sort(gene_list_filtered, decreasing = TRUE)

# GSEA
gsea_go <- gseGO(
  geneList = gene_list_filtered,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Visualize Results
# a. Dotplot of top enriched terms
dotplot(gsea_go, showCategory = 20, title = "GSEA: GO Biological Processes")

# b. Enrichment map plot
gsea_go_enrich <- pairwise_termsim(gsea_go)
emapplot(gsea_go_enrich, showCategory = 10, layout = "kk")

# c. Ridgeplot
ridgeplot(gsea_go) + labs(title = "Distribution of Enrichment Scores")

# d. GSEA plot for a term
# To visualize one enriched pathway in detail:
gseaplot2(gsea_go, geneSetID = "GO:0002449", title = "Lymphocyte Mediated Immunity")



# ----------------------------------------------------------
# Gene Set Enrichment Analysis — Specifically for immune-related gene sets

# MSigDB C7: Immunologic Signatures
msigdb_immune <- msigdbr(species = "Homo sapiens", category = "C7")

# Format TERM2GENE list
msigdb_list <- msigdb_immune[, c("gs_name", "entrez_gene")]

gsea_immune <- GSEA(
  geneList = gene_list_filtered,       # named vector (ENTREZ IDs as names)
  TERM2GENE = msigdb_list,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Dotplot of top pathways
dotplot(gsea_immune, showCategory = 20) + ggtitle("GSEA: Immune-Related Signatures")

# Ridge plot
ridgeplot(gsea_immune) + ggtitle("Enrichment Score Distribution")

# Enrichment map (if enough pathways)
gsea_immune <- pairwise_termsim(gsea_immune)
emapplot(gsea_immune, showCategory = 20, layout = "kk")


# Trying: category = "H" (Hallmark gene sets, immune-related sets like INTERFERON_GAMMA_RESPONSE)
msigdb_immune <- msigdbr(species = "Homo sapiens", category = "H")  

# Format TERM2GENE list
msigdb_list <- msigdb_immune[, c("gs_name", "entrez_gene")]

gsea_immune <- GSEA(
  geneList = gene_list_filtered,       # named vector (ENTREZ IDs as names)
  TERM2GENE = msigdb_list,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Dotplot of top pathways
dotplot(gsea_immune, showCategory = 20) + ggtitle("GSEA: Immune-Related Signatures")

# Ridge plot
ridgeplot(gsea_immune) + ggtitle("Enrichment Score Distribution")

# Enrichment map (if enough pathways)
gsea_immune <- pairwise_termsim(gsea_immune)
emapplot(gsea_immune, showCategory = 20, layout = "kk")
