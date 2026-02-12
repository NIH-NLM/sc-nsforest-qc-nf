/**
 * Silhouette Embedding Dotplot
 *
 * Creates embedding plot (UMAP, t-SNE, etc.) colored by silhouette scores.
 * Helps identify problematic regions in the embedding space.
 */
process viz_dotplot_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    shell '/bin/sh'
    
    input:
    tuple val(meta), path(h5ad)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_dotplot_${meta.embedding}.{html,svg}"),
          emit: plots
    
    script:
    """
    docker run ghcr.io/nih-nlm/scsilhouette:1.0 \
    viz-dotplot \
        --h5ad-path ${h5ad} \
        --embedding-key ${meta.embedding} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
