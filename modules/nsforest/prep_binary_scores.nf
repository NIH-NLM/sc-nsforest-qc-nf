/**
 * Prep Binary Scores Module (Parallelized by Cluster)
 *
 * Computes binary expression scores for each cluster vs all others.
 * Uses ns.pp.prep_binary_scores() on filtered adata from prep_medians.
 */
process prep_binary_scores_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}_${cluster}"
    label 'nsforest'
    
    input:
    tuple val(meta), path(adata_prep), val(cluster)
    
    output:
    tuple val(meta), 
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_binary_scores_*.csv"),
          emit: partial
    
    script:
    """
    nsforest-cli prep-binary-scores \
        --h5ad-path ${adata_prep} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        --cluster-list "${cluster}"
    """
}
