/**
 * Plot Histograms Module
 *
 * Creates histograms of non-zero median and binary score distributions.
 */
process plot_histograms_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(medians_csv), path(binary_scores_csv)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/hist_nonzero_*.svg"),
          emit: histograms
    
    script:
    """
    nsforest-cli plot-histograms \
        --medians-csv ${medians_csv} \
        --binary-scores-csv ${binary_scores_csv} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
