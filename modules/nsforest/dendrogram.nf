/**
 * Dendrogram Module
 *
 * Computes a hierarchical dendrogram over cluster medians and writes the
 * cluster traversal order used to scatter run_nsforest_process.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with organ, first_author, journal, year, author_cell_type, embedding, vid
 *   - h5ad: Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit stats:   tuple(meta, h5ad, cluster_order.csv)  — drives scatter
 * @emit results: tuple(meta, [dendrogram SVG + cluster_order CSV])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{vid}_{cluster_header_safe}_*.{csv,svg}
 */
process dendrogram_process {
    tag "dendrogram_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_cluster_order.csv"),  emit: cluster_order
    tuple val(meta), path("*_cluster_sizes.csv"),  emit: cluster_sizes
    tuple val(meta), path("*_summary_normal.csv"), emit: summary
    tuple val(meta), path("*_dendrogram.svg"),     emit: svg

    script:
    """
    nsforest-cli dendrogram \
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
