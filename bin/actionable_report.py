#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Intersect differential expression results with actionable gene list."
    )
    parser.add_argument("--differential", required=True, help="CSV of DE results")
    parser.add_argument("--actionable", required=True, help="CSV of actionable genes")
    parser.add_argument("--out", required=True, help="Output actionable CSV path")
    parser.add_argument(
        "--summary",
        required=True,
        help="Summary JSON path with key metrics for the run",
    )
    parser.add_argument(
        "--p_adj_cutoff", type=float, default=0.05, help="Adjusted p-value cutoff"
    )
    parser.add_argument(
        "--log2_fc_cutoff", type=float, default=1.0, help="Absolute log2FC cutoff"
    )
    args = parser.parse_args()

    diff = pd.read_csv(args.differential)
    actionable = pd.read_csv(args.actionable)

    if "gene_id" not in diff.columns:
        raise ValueError("Differential results must include gene_id.")
    if "gene_id" not in actionable.columns:
        raise ValueError("Actionable list must include gene_id.")

    filtered = diff[
        (diff["p_adj"] <= args.p_adj_cutoff)
        & (diff["log2_fc"].abs() >= args.log2_fc_cutoff)
    ]
    merged = filtered.merge(actionable, on="gene_id", how="inner")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, index=False)

    summary = {
        "total_genes_tested": len(diff),
        "significant_genes": len(filtered),
        "actionable_hits": len(merged),
        "p_adj_cutoff": args.p_adj_cutoff,
        "log2_fc_cutoff": args.log2_fc_cutoff,
    }
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
