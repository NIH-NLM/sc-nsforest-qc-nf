/**
 * Viz Summary Module
 *
 * Generates an interactive silhouette F-score summary plot combining
 * silhouette scores with NSForest F-scores per cluster.
 * Also writes a dataset-level summary CSV.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:             Map with organ, first_author, journal, year, author_cell_type, embedding, doi, etc.
 *   - silhouette_scores:{prefix}_silhouette_scores.csv
 *   - cluster_summary:  {prefix}_cluster_summary.csv
 *   - annotation:       {prefix}_annotation.json
 *   - nsforest_results: {prefix}_results.csv (or NO_FILE sentinel)
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [summary SVG/HTML/CSV, dataset_summary CSV])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{cluster_header_safe}_*.{csv,svg,html,json}
 */
process viz_summary_process {
    tag "viz_summary_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
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
          path("*.{csv,svg,html,json,log}", optional: true),
          emit: plots

    script:
    """
    scsilhouette viz-summary \
        --silhouette-score-path ${silhouette_scores} \
	--silhouette-score-col ${silhouette_score_col} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
	--dataset-version-id "${meta.dataset_version_id}" \
        --embedding-key "${meta.embedding}" \
	--fscore-path "${nsforest_results}" 
	--doi "${meta.doi}" \
	--collection-name "${meta.collection_name}" \
	--dataset-title "${meta.dataset_title}" \
	--journal "${meta.journal}" \
	--collection_url "${meta.collection_url}" \
	--explorer_url "${meta.explorer_url}" \
	--h5ad_url "${meta.h5ad_url}" 
    """
}
