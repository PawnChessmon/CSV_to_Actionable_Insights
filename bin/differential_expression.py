#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def benjamini_hochberg(p_values):
    """Return FDR-adjusted p-values using Benjamini-Hochberg."""
    p_values = np.array(p_values, dtype=float)
    n = len(p_values)
    order = np.argsort(p_values)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    adjusted = p_values * n / ranks
    adjusted = np.minimum.accumulate(adjusted[order][::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n)
    result[order] = adjusted
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run a simple differential expression test (Welch t-test)."
    )
    parser.add_argument("--counts", required=True, help="Normalized counts CSV")
    parser.add_argument(
        "--metadata", required=True, help="Metadata CSV with sample_id and condition"
    )
    parser.add_argument("--out", required=True, help="Output CSV path")
    args = parser.parse_args()

    counts = pd.read_csv(args.counts)
    metadata = pd.read_csv(args.metadata)

    if "gene_id" not in counts.columns:
        raise ValueError("Counts file must contain a 'gene_id' column.")
    if not {"sample_id", "condition"}.issubset(metadata.columns):
        raise ValueError("Metadata must contain 'sample_id' and 'condition' columns.")

    conditions = metadata["condition"].unique()
    if len(conditions) != 2:
        raise ValueError("Metadata must describe exactly two conditions.")

    cond_a, cond_b = conditions
    samples_a = metadata.loc[metadata["condition"] == cond_a, "sample_id"].tolist()
    samples_b = metadata.loc[metadata["condition"] == cond_b, "sample_id"].tolist()

    missing_a = [s for s in samples_a if s not in counts.columns]
    missing_b = [s for s in samples_b if s not in counts.columns]
    if missing_a or missing_b:
        missing = missing_a + missing_b
        raise ValueError(f"Samples missing from counts matrix: {', '.join(missing)}")

    matrix = counts.set_index("gene_id")
    group_a = matrix[samples_a]
    group_b = matrix[samples_b]

    records = []
    for gene, values_a in group_a.iterrows():
        values_b = group_b.loc[gene]
        a_vals = values_a.values.astype(float)
        b_vals = values_b.values.astype(float)

        test = stats.ttest_ind(a_vals, b_vals, equal_var=False)
        mean_a = a_vals.mean()
        mean_b = b_vals.mean()
        # log2 fold change defined as cond_a - cond_b to match normalization space
        log2_fc = mean_a - mean_b

        records.append(
            {
                "gene_id": gene,
                f"{cond_a}_mean": mean_a,
                f"{cond_b}_mean": mean_b,
                "log2_fc": log2_fc,
                "p_value": float(test.pvalue),
            }
        )

    df = pd.DataFrame.from_records(records)
    df["p_adj"] = benjamini_hochberg(df["p_value"].values)
    df.sort_values("p_value", inplace=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
