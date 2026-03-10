/**
 * Run NSForest Module (Parallelised by Cluster Batch)
 *
 * Runs the NSForest algorithm for a batch of clusters. Each invocation
 * receives the full filtered h5ad plus pre-computed medians and binary
 * scores CSVs, which are loaded into adata.varm before calling NSForest.
 * Scattered over cluster batches so all clusters run in parallel.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map with organ, first_author, year, author_cell_type
 *   - h5ad:        Path to adata_filtered.h5ad
 *   - medians_csv: {prefix}_medians.csv
 *   - binary_csv:  {prefix}_binary_scores.csv
 *   - cluster:     Comma-separated list of cluster names for this batch
 *
 * Output:
 * -------
 * @emit partial: tuple(meta, [partial_results.csv])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_partial_*.csv
 */
process run_nsforest_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}_batch"
    label 'nsforest'

    input:
    tuple val(meta), path(h5ad), path(medians_csv), path(binary_csv), val(cluster)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_*_partial_*.csv", optional: true),
          emit: partial

    script:
    """
    nsforest-cli run-nsforest \
        --h5ad-path ${h5ad} \
        --medians-csv ${medians_csv} \
        --binary-scores-csv ${binary_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --cluster-list '${cluster}' \
        --n-trees 1000 \
        --n-genes-eval 6
    """
}
