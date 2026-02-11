/**
 * Merge Medians Module
 *
 * Combines partial median files from parallel prep_medians jobs.
 * Saves csv + pkl as in DEMO Section 3.
 */
process merge_medians_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(adata_preps), path(partial_csvs)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_medians.{csv,pkl}"),
          emit: complete
    
    script:
    """
    # Take first adata_prep (all identical after prep_medians)
    cp ${adata_preps[0]} outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad

    nsforest-cli merge-medians \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
