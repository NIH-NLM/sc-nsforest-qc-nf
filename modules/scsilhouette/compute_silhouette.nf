/**
 * Compute Silhouette Module
 *
 * Computes per-cell silhouette scores for each cluster using the specified
 * embedding. Saves per-cell scores, per-cluster summary statistics, and
 * an annotation JSON for downstream viz processes.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with organ, first_author, year, author_cell_type, embedding, disease
 *   - h5ad: Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit results: tuple(meta, [silhouette_scores.csv, cluster_summary.csv, annotation.json])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_silhouette_scores.csv
 *                   {organ}_{first_author}_{year}_{cluster_header_safe}_cluster_summary.csv
 *                   {organ}_{first_author}_{year}_{cluster_header_safe}_annotation.json
 */
process compute_silhouette_process {
    tag "compute_silhouette_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("*.{csv,json,log}", optional: true),
          emit: results

    script:
    """
    scsilhouette compute-silhouette \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --embedding-key "${meta.embedding}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
	--dataset-version-id "${meta.dataset_version_id}" \
        --disease "${meta.disease}" \
        --save-scores \
        --save-cluster-summary \
        --save-annotation
    """
}
