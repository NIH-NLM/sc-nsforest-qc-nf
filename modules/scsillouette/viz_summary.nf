/**
 * Silhouette Summary Visualization
 *
 * Creates summary plots of silhouette scores with optional NSForest F-score integration.
 * Generates side-by-side boxplots for comprehensive QC.
 */
process viz_summary_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), 
          path(silhouette_scores), 
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_silhouette_fscore_summary.{html,svg}"),
          emit: plots
    
    script:
    def fscore_flag = nsforest_results.name != 'NO_FILE' ? "--fscore-path ${nsforest_results}" : ""
    """
    scsilhouette viz-summary \
        --silhouette-score-path ${silhouette_scores} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        ${fscore_flag}
    """
}
