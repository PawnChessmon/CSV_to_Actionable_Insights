# CSV to Actionable Insights (demo)

A lightweight Nextflow pipeline showing how to transform a CSV-based cancer expression matrix into differential expression calls and actionable hits.

## What it does
- Normalizes raw counts with log2(CPM + 1) and maps gene IDs to symbols when an annotation table is provided.
- Runs a simple Welch t-test between two conditions, reporting both p-value and BH-adjusted p-value.
- Intersects significant genes with an actionable list.
- Generates PCA, heatmap, volcano, and MA plots (up in red, down in blue, top 5 labeled by p-value).

## Layout
- `main.nf` – pipeline definition.
- `nextflow.config` – defaults and optional conda profile.
- `bin/` – helper scripts executed by processes.
- `data/raw/` – example gene counts + sample metadata.
- `data/reference/` – toy actionable gene list.
- `results/` – created on run; holds outputs.
- `results/plots/` – PCA, heatmap, volcano, and MA plots.

## Quickstart
1. Install Nextflow (https://www.nextflow.io/) and, optionally, `conda` or `mamba`.
2. Run with the bundled toy data:
   ```
   nextflow run main.nf -profile local
   ```
   Or create an isolated environment:
   ```
   nextflow run main.nf -profile conda
   ```
3. Outputs land in `results/`:
   - `preprocessed/normalized_counts.csv`
   - `differential_expression.csv`
   - `actionable_hits.csv`
   - `summary.json`
   - `plots/`:
     - `pca_samples.png`
     - `heatmap_top_genes.png`
     - `volcano.png`
     - `ma_plot.png`

## Customizing inputs
Pass your own files via params:
```
nextflow run main.nf \
  --counts path/to/gene_counts.csv \
  --metadata path/to/sample_metadata.csv \
  --actionable path/to/actionable_genes.csv \
  --annotations path/to/gene_annotations.tsv \
  --outdir my_results
```

Expected formats:
- Counts: `gene_id` column (Ensembl or symbols) plus one column per sample.
- Metadata: `sample_id,condition` with exactly two conditions (e.g., Tumor/Normal).
- Actionable list: `gene_id` plus any annotation columns you like.
- Optional annotations: table with `gene_id` and `gene_symbol` (common headers auto-detected).

## Notes
- Differential expression uses log2 fold-change in the normalized space; adjust thresholds in `bin/actionable_report.py` or `bin/plot_reports.py` if desired.
- Plots use raw p-value <= 0.05 and |log2FC| >= 1 to color up/down (red/blue) and label the top 5 hits by p-value.
- The example data are minimal and meant for workflow illustration, not biological interpretation.
- Omit `results*/`, `work/`, and `.nextflow*` when committing or packaging.
