/**
 * Run NSForest Module (Parallelized by Cluster)
 *
 * Runs NSForest algorithm to identify marker gene combinations.
 * Each cluster processed independently (one vs all).
 */
process run_nsforest_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}_${cluster}"
    label 'nsforest'
    
    input:
    tuple val(meta), path(adata_prep), val(cluster)
    
    output:
    tuple val(meta), 
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_results_*.csv"),
          emit: partial
    
    script:
    """
    nsforest-cli run-nsforest \
        --h5ad-path ${adata_prep} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        --cluster-list "${cluster}" \
        --n-trees 1000 \
        --n-genes-eval 6
    """
}
