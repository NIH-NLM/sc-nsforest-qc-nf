/**
 * Binary Positive Genes Module
 *
 * Lists the genes with a binary score above 0 in each cluster, with the score
 * (one row per cluster and gene), from the binary scores that prep computed.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:       Map with organ, first_author, journal, year, author_cell_type, embedding, dataset_version_id
 *   - binary_csv: binary_scores_ensg_{prefix}.csv (gene-by-cluster)
 *   - h5ad:       the filtered h5ad (for the gene symbols)
 *
 * Output:
 * -------
 * @emit csv: tuple(meta, binary_positive_genes_{prefix}.csv)
 *   Columns: clusterName, gene_ensg, gene_symbol, binary_score
 */
process binary_positive_genes_process {
    tag "binary_positive_genes_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(binary_csv), path(h5ad)

    output:
    tuple val(meta), path("binary_positive_genes_*.csv"), emit: csv

    script:
    """
    nsforest-cli binary-positive-genes \\
        --binary-scores-csv ${binary_csv} \\
        --filtered-h5ad ${h5ad} \\
        --cluster-header "${meta.author_cell_type}" \\
        --organ "${meta.organ}" \\
        --first-author "${meta.first_author}" \\
        --journal "${meta.journal}" \\
        --year "${meta.year}" \\
        --embedding "${meta.embedding}" \\
        --dataset-version-id "${meta.dataset_version_id}"
    """
}
