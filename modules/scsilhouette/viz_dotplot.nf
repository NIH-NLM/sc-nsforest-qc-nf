/**
 * Viz Dotplot Module
 *
 * Generates an embedding scatter plot (UMAP/t-SNE/etc.) coloured by
 * cluster identity, saved as both HTML (interactive) and SVG.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with organ, first_author, year, author_cell_type, embedding
 *   - h5ad: Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [dotplot HTML and SVG])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{cluster_header_safe}_dotplot_{embedding_key}.{html,svg}
 */
process viz_dotplot_process {
    tag "viz_dotplot_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("*.{html,svg,log}", optional: true),
          emit: plots

    script:
    """
    scsilhouette viz-dotplot \
        --h5ad-path ${h5ad} \
        --embedding-key "${meta.embedding}" \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
