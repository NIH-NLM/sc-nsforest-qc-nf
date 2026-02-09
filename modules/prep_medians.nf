/**
 * NSForest Prepare Medians Module
 *
 * Computes median gene expression for individual clusters. Nextflow automatically
 * parallelizes this process - creating one job per cluster from cluster_stats output.
 *
 * Corresponds to DEMO_NS-forest_workflow.ipynb: Section 3
 *
 * Nextflow Parallelization Strategy:
 * ----------------------------------
 * This process demonstrates Nextflow's dataflow parallelization model:
 *
 * 1. SCATTER (Automatic):
 *    - cluster_stats outputs: CSV with N clusters
 *    - Nextflow reads CSV and creates N parallel channels
 *    - Each channel contains: [organ, author, year, cluster_header, h5ad, cluster_name]
 *    - prep_medians_process runs N times in parallel (one per cluster)
 *
 * 2. PARALLEL EXECUTION (Automatic):
 *    - All N instances run concurrently (resource permitting)
 *    - Each computes median for its assigned cluster only
 *    - No coordination needed - pure dataflow
 *
 * 3. GATHER (Automatic):
 *    - Nextflow collects all N partial CSV outputs
 *    - groupTuple() groups by dataset (organ, author, year, cluster_header)
 *    - collectFile() merges partial CSVs into complete median matrix
 *
 * No Python or Bash merge scripts needed - Nextflow handles everything through
 * its dataflow operators.
 *
 * Input:
 * ------
 * @param tuple Cluster-specific job information:
 *   - organ: Organ/tissue
 *   - first_author: First author surname
 *   - year: Publication year
 *   - author_cell_type: Cluster column name
 *   - h5ad: Path to h5ad file
 *   - cluster: Single cluster name to process
 *
 * Output:
 * -------
 * @return tuple: Partial median matrix for this cluster
 *   - organ: Organ/tissue
 *   - first_author: First author surname
 *   - year: Publication year
 *   - author_cell_type: Cluster column name
 *   - partial_csv: Partial median matrix (one row = this cluster)
 *
 * Output Files:
 * -------------
 * {cluster_header}_medians_partial.csv:
 *   - Rows: Single cluster
 *   - Columns: All genes (Ensembl IDs)
 *   - Values: Median expression values
 *
 * These partial files are automatically merged in the workflow to create
 * the complete median matrix with all clusters.
 */
process prep_medians_process {
    tag "${organ}_${first_author}_${year}_${cluster}"
    label 'nsforest'
    
    input:
    tuple val(organ), val(first_author), val(year), val(author_cell_type), 
          path(h5ad), val(cluster)
    
    output:
    tuple val(organ), val(first_author), val(year), val(author_cell_type),
          path("outputs_${organ}_${first_author}_${year}/${author_cell_type}_medians_partial.csv"),
          emit: partial
    
    script:
    """
    nsforest-cli prep-medians \
        --h5ad-path ${h5ad} \
        --cluster-header ${author_cell_type} \
        --organ ${organ} \
        --first-author ${first_author} \
        --year ${year} \
        --cluster-list ${cluster}
    """
}
