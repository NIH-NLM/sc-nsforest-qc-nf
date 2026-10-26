/**
 * Combinatorial Plots Module
 *
 * The scaled dotplot, stacked violin and matrix plot of the global combinatorial NS-Forest results of the clusters on all
 * the cells, as in the Zenodo deposit (dotplot_global_cluster_combinatorial_on_global.png): one row per cluster, and under
 * each cluster its best combination of subclass and cluster markers; each as PNG, SVG and HTML. The rows follow the order
 * of the dendrogram file (Jorstad Fig 5B).
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map of the published cluster header
 *   - h5ad:        the h5ad with the cluster column
 *   - results_csv: the NS-Forest results of the paper for the clusters (the same file as plots_process)
 *   - dend_order:  CSV with the cluster order of the paper (column CrossArea_cluster_order)
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [dotplot, stacked violin and matrix plot as PNG, SVG and HTML])
 */
process plots_combinatorial_process {
    tag "plots_combinatorial_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), path(results_csv), path(dend_order, stageAs: 'in_dend.csv')

    output:
    tuple val(meta),
          path("*.{svg,html,png,log}", optional: true),
          emit: plots

    script:
    """
    nsforest-cli plots-combinatorial \
        --h5ad-path ${h5ad} \
        --results-csv ${results_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}" \
        --dend-order-csv ${dend_order}
    """
}
