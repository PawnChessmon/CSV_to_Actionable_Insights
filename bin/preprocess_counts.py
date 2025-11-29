#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Return log2(CPM + 1) transformed counts."""
    counts_per_million = counts.div(counts.sum()) * 1_000_000
    return (counts_per_million + 1).applymap(lambda x: np.log2(x))


def main():
    parser = argparse.ArgumentParser(
        description="Normalize gene counts and align with sample metadata."
    )
    parser.add_argument("--counts", required=True, help="CSV with gene counts")
    parser.add_argument("--metadata", required=True, help="CSV with sample metadata")
    parser.add_argument(
        "--annotations",
        required=False,
        help="Optional gene annotation table with gene ID and gene symbol",
    )
    parser.add_argument("--out", required=True, help="Output CSV path")
    args = parser.parse_args()

    counts_path = Path(args.counts)
    metadata_path = Path(args.metadata)

    counts = pd.read_csv(counts_path)
    metadata = pd.read_csv(metadata_path)

    if "gene_id" not in counts.columns:
        raise ValueError("Counts file must contain a 'gene_id' column.")
    if not {"sample_id", "condition"}.issubset(metadata.columns):
        raise ValueError("Metadata must contain 'sample_id' and 'condition' columns.")

    if args.annotations:
        ann = pd.read_csv(
            args.annotations, sep=None, engine="python"
        )  # auto-detect comma/tsv
        # heuristics for id and symbol columns
        id_cols = [c for c in ann.columns if c.lower().replace(" ", "_") in ["gene_id", "geneid", "id"]]
        sym_cols = [
            c
            for c in ann.columns
            if c.lower().replace(" ", "_")
            in ["gene_symbol", "symbol", "associated_gene_name", "gene_name"]
        ]
        if id_cols and sym_cols:
            ann = ann[[id_cols[0], sym_cols[0]]].rename(
                columns={id_cols[0]: "gene_id", sym_cols[0]: "gene_symbol"}
            )
            ann = ann.dropna(subset=["gene_symbol"])
            ann_map = ann.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]
            counts["gene_symbol"] = counts["gene_id"].map(ann_map)
            # prefer symbol; fallback to original id when missing
            counts["gene_id"] = counts["gene_symbol"].fillna(counts["gene_id"])
            counts = counts.drop(columns=["gene_symbol"])
            # remove duplicate symbols to keep unique index later
            counts = counts.drop_duplicates(subset=["gene_id"])

    # Keep only samples present in metadata and maintain the same ordering
    sample_ids = metadata["sample_id"].tolist()
    missing = [s for s in sample_ids if s not in counts.columns]
    if missing:
        raise ValueError(f"Samples missing from counts matrix: {', '.join(missing)}")

    counts_matrix = counts.set_index("gene_id")[sample_ids]
    normalized = log2_cpm(counts_matrix)

    normalized.insert(0, "gene_id", normalized.index)
    normalized.reset_index(drop=True, inplace=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
