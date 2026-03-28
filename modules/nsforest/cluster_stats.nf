/**
 * Cluster Statistics Module
 *
 * Computes basic cluster statistics (cell counts, percentages).
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with organ, first_author, year, author_cell_type
 *   - h5ad: Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit results: tuple(meta, [cluster_statistics CSV])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{cluster_header}_{vid}_cluster_statistics.csv
 */
process cluster_stats_process {
    tag "cluster_stats_${meta.organ}_${meta.first_author}_${meta.year}_${meta.journal}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("*.{csv,svg,html,log}", optional: true),
          emit: results

    script:
    """
    nsforest-cli cluster-stats \
        --h5ad-path "${h5ad}" \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
	--embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
