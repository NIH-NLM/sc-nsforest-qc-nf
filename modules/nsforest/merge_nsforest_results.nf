/**
 * Merge NSForest Results Module
 *
 * Combines partial NSForest results from parallel run_nsforest jobs.
 * Saves csv + pkl.
 */
process merge_nsforest_results_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(partial_csvs)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_results.{csv,pkl}"),
          emit: complete
    
    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
