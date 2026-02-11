/**
 * Filter AnnData Module
 *
 * Filters by disease status, tissue, and minimum cluster size.
 * Creates before/after dendrograms and statistics for comparison.
 */
process filter_adata_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(h5ad)
    
    output:
    tuple val(meta), 
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_filtered.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}*.{csv,svg}"),
          emit: results
    
    script:
    def filter_normal_flag = meta.filter == "True" ? "--filter-normal" : ""
    """
    nsforest-cli filter-adata \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        ${filter_normal_flag} \
        --tissue ${meta.tissue} \
        --min-cluster-size 5
    """
}
