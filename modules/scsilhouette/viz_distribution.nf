/**
 * Viz Distribution Module
 *
 * Generates distribution plots of cluster cell counts (raw and log10)
 * overlaid with mean/median silhouette scores per cluster.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:            Map with organ, first_author, year, author_cell_type
 *   - silhouette_scores: {prefix}_silhouette_scores.csv
 *   - cluster_summary: {prefix}_cluster_summary.csv
 *   - annotation:      {prefix}_annotation.json
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [distribution HTML and SVG])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_distribution_*.{html,svg}
 */
process viz_distribution_process {
    tag "viz_distribution_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${dataset_version_id}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta),
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation)

    output:
    tuple val(meta),
          path("*.{csv,svg,html,json}", optional: true),
          emit: plots

    script:
    """
    scsilhouette viz-distribution \
        --cluster-summary-path ${cluster_summary} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
	--embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
