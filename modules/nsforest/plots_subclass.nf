/**
 * Subclass Plots Module
 *
 * The stacked violin and the dotplot of the NS-Forest results of the subclasses on all the cells, as in the
 * Zenodo deposit (stacked_violin_global_subclass.png, dotplot_global_subclass.png): one row per subclass, and
 * under each subclass its own markers. The rows follow the order of the dendrogram file (Jorstad Fig 5B).
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map of the subclass header
 *   - h5ad:        the h5ad with the subclass and the cluster columns
 *   - results_csv: NS-Forest results of the subclasses
 *   - dend_order:  CSV with the cluster order of the paper (column CrossArea_cluster_order)
 * @param child_header: obs column of the clusters named in the dendrogram file
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [stacked violin and dotplot as SVG])
 */
process plots_subclass_process {
    tag "plots_subclass_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), path(results_csv, stageAs: 'in_results.csv'), path(dend_order, stageAs: 'in_dend.csv')
    val child_header   // the cluster header of the dendrogram file (a module does not see the params of main.nf)

    output:
    tuple val(meta),
          path("*.{svg,html,log}", optional: true),
          emit: plots

    script:
    """
    nsforest-cli plots-subclass \
        --h5ad-path ${h5ad} \
        --results-csv ${results_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --child-header "${child_header}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}" \
        --dend-order-csv ${dend_order}
    """
}
