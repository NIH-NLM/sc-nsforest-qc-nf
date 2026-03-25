/**
 * Compute Summary Statistics
 *
 * Creates dataset-level summary statistics from cluster summaries.
 * Computes median-of-medians and other aggregate metrics across all clusters.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:             Map with organ, first_author, year, author_cell_type, embedding, doi, etc.
 *   - silhouette_scores:{prefix}_silhouette_scores.csv
 *   - cluster_summary:  {prefix}_cluster_summary.csv
 *   - annotation:       {prefix}_annotation.json
 *   - nsforest_results: {prefix}_results.csv (or NO_FILE sentinel)
 *
 * Output:
 * -------
 * @emit summary: tuple(meta, {prefix}_dataset_summary.csv)
 *   Contains: organ, first_author, year, cluster_header, n_clusters, n_cells,
 *             median/mean/std silhouette, quality tier counts, median/mean F-score,
 *             doi, collection_name, dataset_title, journal
 */
process compute_summary_stats_process {
    tag "compute_summary_stats_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta),
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)

    output:
    tuple val(meta),
          path("*.csv"),
          emit: summary

    script:
    """
    scsilhouette compute-summary-stats \
        --cluster-summary "${cluster_summary}" \
        --nsforest-results "${nsforest_results}" \
	--metadata "${annotation}" \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --embedding-key "${meta.embedding}" 
    """
}
