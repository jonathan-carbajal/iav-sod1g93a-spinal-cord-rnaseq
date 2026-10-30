# ADVANCED PREVIEW SCRIPT -- NOT the HOSA submission.
# This demonstrates techniques beyond the current project: full QC plots,
# a volcano plot, a top-gene heatmap, and gene set enrichment analysis (GSEA).
# Every number here is computed for real from the actual data -- nothing is
# invented -- but this file is a "here's where this could go" exploration,
# not something to present as independently-performed current work.

library(limma)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/bioresearch-workspace/iav-sod1g93a-spinal-cord-rnaseq")

## ---- Rebuild the canonical pipeline (same as iav_sod1g93a_analysis.R) ----

dat <- read.csv("RNAseq_Mouse_Gene.csv")
rownames(dat) <- dat$X
dat$X <- NULL
expr <- as.matrix(dat)

metadata <- data.frame(
  sample = colnames(expr),
  condition = c("Uninfected", "Uninfected", "Uninfected",
                "Infected", "Infected", "Infected")
)
rownames(metadata) <- colnames(expr)
metadata$condition <- factor(metadata$condition, levels = c("Uninfected", "Infected"))

design <- model.matrix(~ condition, data = metadata)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "conditionInfected", number = Inf, sort.by = "P")

cat("Genes with adj.P.Val < 0.05:", sum(results$adj.P.Val < 0.05), "of", nrow(results), "\n")

## ---- 1. QC plots (guide section 10, never finished in the main project) ----

png("advanced_preview/figures/01_boxplot.png", width = 1400, height = 900, res = 150)
boxplot(expr, las = 2, main = "Expression distributions by sample",
        ylab = "Normalized log-expression", cex.axis = 0.7)
dev.off()

pca <- prcomp(t(expr), scale. = FALSE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                      condition = metadata$condition, sample = rownames(metadata))
pct_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA of samples",
       x = paste0("PC1 (", pct_var[1], "% variance)"),
       y = paste0("PC2 (", pct_var[2], "% variance)")) +
  theme_minimal(base_size = 13)
ggsave("advanced_preview/figures/02_pca.png", p_pca, width = 7, height = 5, dpi = 150)

png("advanced_preview/figures/03_sample_correlation_heatmap.png",
    width = 1000, height = 900, res = 150)
pheatmap(cor(expr, method = "pearson"),
         main = "Sample-to-sample correlation",
         annotation_col = data.frame(condition = metadata$condition,
                                     row.names = rownames(metadata)))
dev.off()

## ---- 2. Volcano plot ----
## Being careful with labeling here: NO gene in this project reached
## adj.P.Val < 0.05. The "genes of interest" below are flagged only at a
## loose, UNADJUSTED p < 0.01 threshold, purely to show the technique --
## this is explicitly not a claim of statistical significance.

results$flag <- ifelse(results$P.Value < 0.01, "nominal p<0.01 (not adj.-significant)", "not notable")

p_volcano <- ggplot(results, aes(logFC, -log10(P.Value), color = flag)) +
  geom_point(alpha = 0.6, size = 1.3) +
  scale_color_manual(values = c("nominal p<0.01 (not adj.-significant)" = "firebrick",
                                 "not notable" = "grey70")) +
  labs(title = "Volcano plot: Infected vs. Uninfected",
       subtitle = "No gene reached adj.P.Val < 0.05 in this 3-vs-3 comparison",
       x = "log2 fold change (Infected - Uninfected)", y = "-log10(raw P.Value)",
       color = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave("advanced_preview/figures/04_volcano.png", p_volcano, width = 7, height = 5.5, dpi = 150)

## ---- 3. Heatmap of the top-ranked genes (by raw p-value, NOT "significant") ----

top_genes <- rownames(results)[order(results$P.Value)][1:30]
mat <- expr[top_genes, ]
mat_scaled <- t(scale(t(mat)))  # z-score each gene across samples for visualization

png("advanced_preview/figures/05_top30_genes_heatmap.png",
    width = 1100, height = 1300, res = 150)
pheatmap(mat_scaled,
         main = "Top 30 genes by raw p-value (NOT adjusted-significant)",
         annotation_col = data.frame(condition = metadata$condition,
                                     row.names = rownames(metadata)),
         fontsize_row = 7)
dev.off()

## ---- 4. GSEA: pathway-level analysis using the FULL ranked gene list ----
## This is the key "advanced" idea: instead of needing a cutoff list of
## significant genes (which is empty here), GSEA ranks ALL ~14,886 genes by
## their t-statistic and asks whether a whole pathway's (GO Biological
## Process term's) genes are skewed toward one end of that ranking -- a
## coordinated shift no single gene reaches alone.
## Uses GO terms via org.Mm.eg.db (installed locally) instead of MSigDB
## Hallmark sets, because MSigDB's own gene-set download host (Zenodo) is
## unreachable from this environment -- documented here rather than silently
## worked around, same as any other real-world data-access limitation.

ranks <- results$t
names(ranks) <- rownames(results)
ranks <- sort(ranks, decreasing = TRUE)
ranks <- ranks[!duplicated(names(ranks))]

set.seed(1)
gsea_res <- gseGO(geneList = ranks, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                   ont = "BP", minGSSize = 15, maxGSSize = 500,
                   pvalueCutoff = 1, eps = 0, seed = TRUE)
gsea_df <- as.data.frame(gsea_res)
gsea_df <- gsea_df[order(gsea_df$pvalue), ]

write.csv(gsea_df[, c("ID", "Description", "pvalue", "p.adjust", "NES", "setSize")],
          "advanced_preview/gsea_go_bp_results.csv", row.names = FALSE)

cat("\nTop 10 GO Biological Process terms by raw p-value (GSEA):\n")
print(gsea_df[1:10, c("Description", "pvalue", "p.adjust", "NES")])
cat("\nGO terms with p.adjust < 0.05:", sum(gsea_df$p.adjust < 0.05, na.rm = TRUE),
    "of", nrow(gsea_df), "\n")

top20 <- gsea_df[1:20, ]
top20$Description <- factor(top20$Description, levels = rev(top20$Description))

p_gsea <- ggplot(top20, aes(NES, Description, fill = pvalue < 0.05)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "grey70"),
                     labels = c(`TRUE` = "raw p < 0.05", `FALSE` = "not nominally sig."),
                     name = NULL) +
  labs(title = "Top 20 GO Biological Process terms by GSEA",
       subtitle = "NES = normalized enrichment score; see p.adjust in the CSV before over-interpreting",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_minimal(base_size = 10)
ggsave("advanced_preview/figures/06_gsea_top_pathways.png", p_gsea, width = 10, height = 7, dpi = 150)

cat("\nDone. Figures in advanced_preview/figures/, table in advanced_preview/gsea_go_bp_results.csv\n")
