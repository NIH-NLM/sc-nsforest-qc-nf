/**
 * Merge NSForest Results Module
 *
 * Collects all partial results CSVs from scattered run_nsforest_process
 * invocations and merges them into a single results file.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:         Map with organ, first_author, year, author_cell_type
 *   - partial_csvs: List of partial results CSVs (one per cluster batch)
 *
 * Output:
 * -------
 * @emit complete: tuple(meta, [results.csv, results.pkl])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_results.{csv,pkl}
 */
process merge_nsforest_results_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "*.{csv,pkl}"

    input:
    tuple val(meta), path(partial_csvs, stageAs: "partial_??.csv")

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_*.{csv,pkl}", optional: true),
          emit: complete

    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
