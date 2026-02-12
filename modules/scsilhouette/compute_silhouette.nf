/**
 * Compute Silhouette Scores
 *
 * Computes silhouette scores for single-cell clusters using embeddings.
 * Outputs per-cell scores, cluster summary statistics, and metadata.
 */
process compute_silhouette_process {
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
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_silhouette_scores.csv"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_cluster_summary.csv"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_annotation.json"),
          emit: results
    
    script:
    """
    docker run ghcr.io/nih-nlm/scsilhouette:1.0 \
        compute-silhouette \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --embedding-key ${meta.embedding} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
