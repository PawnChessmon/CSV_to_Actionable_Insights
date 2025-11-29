nextflow.enable.dsl=2

/*
 * CSV to actionable insights demo using a small cancer differential expression set.
 * Steps:
 * 1) Normalize counts
 * 2) Run simple differential expression
 * 3) Flag actionable genes
 */

params.counts    = params.counts    ?: "data/raw/gene_counts.csv"
params.metadata  = params.metadata  ?: "data/raw/sample_metadata.csv"
params.actionable = params.actionable ?: "data/reference/actionable_genes.csv"
params.annotations = params.annotations ?: "data/NEW_data/annotations.csv"
params.outdir    = params.outdir    ?: "results"

process PREPROCESS_COUNTS {
    tag "$counts_file.baseName"
    publishDir "${params.outdir}/preprocessed", mode: 'copy'

    input:
    path counts_file
    path metadata_file
    path annotations_file

    output:
    path "normalized_counts.csv"

    script:
    def annArg = annotations_file ? "--annotations ${annotations_file}" : ""
    """
    preprocess_counts.py \\
        --counts $counts_file \\
        --metadata $metadata_file \\
        --out normalized_counts.csv \\
        $annArg
    """
}

process DIFFERENTIAL_EXPRESSION {
    tag "$normalized_file.baseName"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path normalized_file
    path metadata_file

    output:
    path "differential_expression.csv"

    script:
    """
    differential_expression.py \\
        --counts $normalized_file \\
        --metadata $metadata_file \\
        --out differential_expression.csv
    """
}

process ACTIONABLE_REPORT {
    tag "$diff_file.baseName"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path diff_file
    path actionable_file

    output:
    path "actionable_hits.csv"
    path "summary.json"

    script:
    """
    actionable_report.py \\
        --differential $diff_file \\
        --actionable $actionable_file \\
        --out actionable_hits.csv \\
        --summary summary.json
    """
}

process PLOT_REPORTS {
    tag "$diff_file.baseName"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path normalized_file
    path metadata_file
    path diff_file

    output:
    path "pca_samples.png"
    path "heatmap_top_genes.png"
    path "volcano.png"
    path "ma_plot.png"

    script:
    """
    plot_reports.py \\
        --counts $normalized_file \\
        --metadata $metadata_file \\
        --differential $diff_file
    """
}

workflow {
    counts_ch = Channel.fromPath(params.counts)
    actionable_ch = Channel.fromPath(params.actionable)
    annotations_ch = Channel.fromPath(params.annotations)

    metadata_for_preprocess = Channel.fromPath(params.metadata)
    metadata_for_diff = Channel.fromPath(params.metadata)
    metadata_for_plots = Channel.fromPath(params.metadata)

    normalized_ch = PREPROCESS_COUNTS(counts_ch, metadata_for_preprocess, annotations_ch)

    differential_ch = DIFFERENTIAL_EXPRESSION(normalized_ch, metadata_for_diff)
    ACTIONABLE_REPORT(differential_ch, actionable_ch)

    PLOT_REPORTS(normalized_ch, metadata_for_plots, differential_ch)
}
