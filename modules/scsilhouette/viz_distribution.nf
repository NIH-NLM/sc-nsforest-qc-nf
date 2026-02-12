/**
 * Silhouette Distribution Plots
 *
 * Creates distribution plots showing cluster sizes vs silhouette scores.
 * Generates both log10 and raw scale versions to identify outliers.
 */
process viz_distribution_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    shell '/bin/sh'
    
    input:
    tuple val(meta), 
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_distribution_*.{html,svg}"),
          emit: plots
    
    script:
    """
    docker run ghcr.io/nih-nlm/scsilhouette:1.0 \
    viz-distribution \
        --cluster-summary-path ${cluster_summary} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
