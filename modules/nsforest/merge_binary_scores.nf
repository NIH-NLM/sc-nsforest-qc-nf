/**
 * Merge Binary Scores Module
 *
 * Combines partial binary scores files from parallel prep_binary_scores jobs.
 * Saves csv + pkl.
 */
process merge_binary_scores_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(adata_preps, stageAs: 'adata_prep_*.h5ad'), path(partial_csvs)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_binary_scores.{csv,pkl}"),
          emit: complete
    
    script:
    """
    # Take first adata_prep (all identical after prep_binary_scores)
    mkdir -p outputs_${meta.organ}_${meta.first_author}_${meta.year}
    cp ${adata_preps[0]} outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad

    nsforest-cli merge-binary-scores \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
