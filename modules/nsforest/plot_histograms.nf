/**
 * Plot Histograms Module
 *
 * Creates histograms of non-zero median and binary score distributions
 * to help assess marker gene signal quality before running NSForest.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:             Map with organ, first_author, year, author_cell_type
 *   - medians_csv:      {prefix}_medians.csv
 *   - binary_scores_csv:{prefix}_binary_scores.csv
 *
 * Output:
 * -------
 * @emit histograms: tuple(meta, [hist_nonzero_*.svg])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_hist_nonzero_*.svg
 */
process plot_histograms_process {
    tag "plot_histograms_${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(medians_csv), path(binary_scores_csv)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}*.{svg,html,log}"),
          emit: histograms

    script:
    """
    nsforest-cli plot-histograms \
        --medians-csv ${medians_csv} \
        --binary-scores-csv ${binary_scores_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
