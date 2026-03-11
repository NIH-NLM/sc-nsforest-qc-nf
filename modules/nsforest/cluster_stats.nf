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
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header}_cluster_statistics.csv
 */
process cluster_stats_process {
    tag "cluster_stats_${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_*.{csv,svg,html,log}", optional: true),
          emit: results

    script:
    """
    nsforest-cli cluster-stats \\
        --h5ad-path ${h5ad} \\
        --cluster-header "${meta.author_cell_type}" \\
        --organ "${meta.organ}" \\
        --first-author "${meta.first_author}" \\
        --year "${meta.year}"
    """
}
