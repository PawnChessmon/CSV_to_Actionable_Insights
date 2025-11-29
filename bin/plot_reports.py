#!/usr/bin/env python3
import argparse
import os
import tempfile
from pathlib import Path

# Ensure matplotlib can write its cache in sandboxed environments
mpl_config = Path(tempfile.mkdtemp(prefix="mplconfig_"))
os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


def make_pca(counts: pd.DataFrame, metadata: pd.DataFrame, out_path: Path):
    """Plot sample-level PCA."""
    # Drop gene_id and transpose to samples x genes
    matrix = counts.set_index("gene_id").T
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(matrix.values)
    pc_df = pd.DataFrame(coords, index=matrix.index, columns=["PC1", "PC2"])
    pc_df = pc_df.join(metadata.set_index("sample_id"))

    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=pc_df,
        x="PC1",
        y="PC2",
        hue="condition",
        s=90,
        edgecolor="black",
    )
    plt.title("PCA of samples")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_heatmap(counts: pd.DataFrame, metadata: pd.DataFrame, out_path: Path, top_n: int):
    """Heatmap of top variable genes."""
    matrix = counts.set_index("gene_id")
    # Variance across samples
    top_genes = matrix.var(axis=1).sort_values(ascending=False).head(top_n).index
    mat_top = matrix.loc[top_genes]
    # Reorder columns by condition
    ordered_samples = (
        metadata.sort_values("condition")["sample_id"].tolist()
        if "sample_id" in metadata.columns
        else matrix.columns.tolist()
    )
    mat_top = mat_top[ordered_samples]

    plt.figure(figsize=(8, 10))
    sns.heatmap(mat_top, cmap="magma", yticklabels=True)
    plt.title("Top variable genes (log2 CPM)")
    plt.xlabel("Samples")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_volcano(diff: pd.DataFrame, out_path: Path):
    """Volcano plot of differential results."""
    df = diff.copy()
    df["neg_log10_p"] = -np.log10(df["p_value"].replace(0, np.nextafter(0, 1)))
    up = (df["p_value"] <= 0.05) & (df["log2_fc"] >= 1)
    down = (df["p_value"] <= 0.05) & (df["log2_fc"] <= -1)

    plt.figure(figsize=(6, 5))
    plt.scatter(df.loc[down, "log2_fc"], df.loc[down, "neg_log10_p"], c="royalblue", s=25, alpha=0.8, label="Down")
    plt.scatter(df.loc[up, "log2_fc"], df.loc[up, "neg_log10_p"], c="crimson", s=25, alpha=0.8, label="Up")
    plt.scatter(df.loc[~(up | down), "log2_fc"], df.loc[~(up | down), "neg_log10_p"], c="lightgray", s=20, alpha=0.6, label="NS")

    top_labels = df.sort_values("p_value").head(5)
    for _, row in top_labels.iterrows():
        plt.text(row["log2_fc"], row["neg_log10_p"], row["gene_id"], fontsize=7, ha="right")

    plt.axvline(1, color="black", ls="--", lw=1)
    plt.axvline(-1, color="black", ls="--", lw=1)
    plt.axhline(-np.log10(0.05), color="black", ls="--", lw=1)
    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 p-value")
    plt.title("Volcano plot")
    plt.legend()
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_ma(diff: pd.DataFrame, counts: pd.DataFrame, out_path: Path):
    """MA plot using normalized counts."""
    matrix = counts.set_index("gene_id")
    avg_expr = matrix.mean(axis=1)

    df = diff.set_index("gene_id").copy()
    df["avg_expr"] = avg_expr

    plt.figure(figsize=(6, 5))
    down = (df["p_value"] <= 0.05) & (df["log2_fc"] <= -1)
    up = (df["p_value"] <= 0.05) & (df["log2_fc"] >= 1)
    plt.scatter(df.loc[down, "avg_expr"], df.loc[down, "log2_fc"], c="royalblue", s=25, alpha=0.8, label="Down")
    plt.scatter(df.loc[up, "avg_expr"], df.loc[up, "log2_fc"], c="crimson", s=25, alpha=0.8, label="Up")
    plt.scatter(df.loc[~(up | down), "avg_expr"], df.loc[~(up | down), "log2_fc"], c="lightgray", s=20, alpha=0.6, label="NS")

    top_labels = df.sort_values("p_value").head(5).reset_index()
    for _, row in top_labels.iterrows():
        plt.text(row["avg_expr"], row["log2_fc"], row["gene_id"], fontsize=7, ha="right")

    plt.axhline(0, color="black", lw=1)
    plt.xlabel("Average expression (log2 CPM)")
    plt.ylabel("log2 fold change")
    plt.title("MA plot")
    plt.legend()
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Create PCA, heatmap, volcano, and MA plots.")
    parser.add_argument("--counts", required=True, help="Normalized counts CSV")
    parser.add_argument("--metadata", required=True, help="Metadata CSV")
    parser.add_argument("--differential", required=True, help="Differential expression CSV")
    parser.add_argument("--top_n_heatmap", type=int, default=30, help="Top variable genes to plot")
    args = parser.parse_args()

    counts = pd.read_csv(args.counts)
    metadata = pd.read_csv(args.metadata)
    diff = pd.read_csv(args.differential)

    make_pca(counts, metadata, Path("pca_samples.png"))
    make_heatmap(counts, metadata, Path("heatmap_top_genes.png"), args.top_n_heatmap)
    make_volcano(diff, Path("volcano.png"))
    make_ma(diff, counts, Path("ma_plot.png"))


if __name__ == "__main__":
    main()
