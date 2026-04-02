/**
 * Viz 2D_Projection Module
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
 * @emit plots: tuple(meta, [2D_projection HTML and SVG])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{cluster_header_safe}_{embedding_key}_2D_projection.{html,svg}
 */
process viz_2D_projection_process {
    tag "viz_2D_projection_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
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
    scsilhouette viz-2D-projection \
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
