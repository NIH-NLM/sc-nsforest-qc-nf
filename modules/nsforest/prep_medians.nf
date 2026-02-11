/**
 * Prep Medians Module (Parallelized by Cluster)
 *
 * Corresponds to DEMO_NS-forest_workflow.py: Section 3
 * Uses ns.pp.prep_medians() to compute median expression per cluster.
 */
process prep_medians_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}_${cluster}"
    label 'nsforest'
    
    input:
    tuple val(meta), path(h5ad), val(cluster)
    
    output:
    tuple val(meta), 
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_medians_*.csv"),
          emit: partial
    
    script:
    """
    nsforest-cli prep-medians \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        --cluster-list "${cluster}"
    """
}
