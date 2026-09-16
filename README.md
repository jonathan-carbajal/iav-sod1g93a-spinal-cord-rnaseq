# iav-sod1g93a-spinal-cord-rnaseq

Secondary analysis of public RNA-seq data, done as a learning project for a HOSA Research Poster.

## Research question

Does influenza A virus (IAV) infection alter gene expression in the lumbar spinal cord of `SOD1G93A` mice (an ALS mouse model), 30 days after infection? Secondary analysis of GEO accession [GSE342036](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342036), comparing 3 uninfected vs. 3 infected mice.

## Contents

- `RNAseq_Mouse_Gene.csv` — expression matrix (14,886 genes × 6 samples).
- `iav_sod1g93a_analysis.R` — analysis script: builds the expression matrix and sample metadata, fits a linear model per gene with `limma` (`lmFit`/`eBayes`), and ranks results with `topTable`.

## Status

Sections 1–9 of the underlying workflow (data loading through differential-expression results) are complete: no individual gene reached the BH-adjusted significance threshold (adj.P.Val < 0.05) in this preliminary comparison — a result consistent with the small sample size (n=3 per group), not evidence that infection has no effect. QC plots (boxplot/PCA/correlation heatmap) and a still-open discrepancy between this file's 14,886 genes and the 12,503 filtered genes reported in the public GEO record are not yet addressed.
