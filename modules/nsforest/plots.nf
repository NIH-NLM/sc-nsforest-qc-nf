/**
 * Plots Module
 *
 * Creates boxplots, scatter plots, and expression dotplots for NSForest
 * marker genes. Maps Ensembl IDs to gene symbols using cell-kn gene mapping.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map with organ, first_author, year, author_cell_type
 *   - h5ad:        Path to adata_filtered.h5ad
 *   - results_csv: {prefix}_results.csv from merge_nsforest_results_process
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [boxplots, scatter plots, dotplots as SVG/HTML])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_*.{svg,html}
 */
process plots_process {
    tag "plots_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), path(results_csv)

    output:
    tuple val(meta),
          path("*.{svg,html,log}", optional: true),
          emit: plots

    script:
    """
    nsforest-cli plots \
        --h5ad-path ${h5ad} \
        --results-csv ${results_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
