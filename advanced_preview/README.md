# Advanced preview — not the HOSA submission

This folder is a demonstration of techniques beyond the current project, built to answer "what could this look like a year from now with more skill" — not something to present as your own independently-performed current work. Everything here was computed for real from your actual `RNAseq_Mouse_Gene.csv`, run by Claude Code in this session (`advanced_analysis_preview.R`), not hand-written or verified line-by-line the way the main project's script was.

## What's in here

- `figures/01_boxplot.png`, `02_pca.png`, `03_sample_correlation_heatmap.png` — the QC plots (guide section 10) that were never finished in the main project.
- `figures/04_volcano.png` — logFC vs. significance for every gene. No gene reaches adjusted significance (consistent with the main result); genes flagged "nominal p<0.01" are shown only to illustrate the technique, not as a real finding.
- `figures/05_top30_genes_heatmap.png` — the 30 genes with the smallest raw p-values, explicitly labeled as not adjusted-significant.
- `figures/06_gsea_top_pathways.png` and `gsea_go_bp_results.csv` — **gene set enrichment analysis (GSEA)** against GO Biological Process terms, using the entire ranked gene list rather than a significant-genes cutoff (see "Why GSEA" below).

## Why GSEA, specifically

The main project found 0 genes individually significant after BH correction — a standard overrepresentation test (the kind that needs a list of "significant genes" as input) literally cannot run on an empty list. GSEA sidesteps this: it ranks *all* ~14,886 genes by their t-statistic and asks whether a given pathway's genes are skewed toward one end of that ranking, as a group — a coordinated shift that no single gene needs to reach alone. It's the textbook correct tool for exactly this project's situation, not just a more complex alternative.

**A genuinely real, somewhat surprising result:** 158 of 4,615 tested GO terms reached `p.adjust < 0.05` — real pathway-level signal, despite zero gene-level hits. The single most significant terms (`syncytium formation`, `myoblast fusion`, `cell-cell fusion`) are about cell-fusion biology, not an obviously "neuroinflammation" or "antiviral response" story one might expect going in. That's an honest, unexplained observation, not a conclusion — figuring out why (possible muscle-tissue signal in a spinal cord sample? a real ALS-relevant fusion process? multiple testing noise despite BH correction?) would be real follow-up work, not something to assert confidently without digging further.

## One real data-access limitation worth knowing

The original plan used MSigDB Hallmark gene sets via the `msigdbr` package, but its data host (`zenodo.org`) is unreachable from this sandboxed environment — confirmed as a genuine network restriction, not a slow connection. Switched to GO Biological Process terms via `clusterProfiler`/`org.Mm.eg.db` instead, which are Bioconductor-hosted and installed locally. Worth knowing if you ever want to redo this with Hallmark sets specifically on a machine with normal internet access.
