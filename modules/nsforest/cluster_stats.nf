/**
 * NSForest Cluster Statistics Module
 *
 * Computes basic statistics for each cluster including cell counts and percentages.
 * The output cluster list drives Nextflow's automatic parallelization of downstream
 * processes (prep_medians, run_nsforest).
 *
 * Corresponds to DEMO_NS-forest_workflow.ipynb: Section 2
 *
 * Description:
 * ------------
 * This process analyzes the distribution of cells across clusters and generates
 * a summary table. This is a critical step because:
 * 
 * 1. Quality Control: Identifies clusters with very few cells (potential artifacts)
 * 2. Parallelization: The cluster list enables Nextflow to scatter subsequent
 *    computationally intensive steps (median calculation, NSForest) across clusters
 * 3. Documentation: Provides clear record of cluster composition
 *
 * Nextflow Dataflow Integration:
 * ------------------------------
 * The output CSV's 'cluster' column is read by the workflow to create parallel
 * channels for per-cluster processing. For a dataset with N clusters:
 * 
 *   cluster_stats → [cluster1, cluster2, ..., clusterN] → N parallel jobs
 *
 * This happens automatically through Nextflow's flatMap and groupTuple operators
 * in the main workflow - no explicit parallel code needed here.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map containing dataset metadata
 *       * organ: Organ/tissue type
 *       * first_author: First author surname
 *       * year: Publication year
 *       * author_cell_type: Cluster column name in AnnData.obs
 *       * embedding: Embedding key (X_umap, X_tsne, etc.)
 *       * disease: Disease state
 *       * filter: Filter normal cells flag
 *       * tissue: Tissue type
 *   - h5ad: Path to h5ad file
 *
 * Output:
 * -------
 * @return tuple:
 *   - meta: Same metadata map (passed through)
 *   - h5ad: Same h5ad file (passed through for downstream processes)
 *   - stats_csv: Cluster statistics CSV file
 *
 * Output Files:
 * -------------
 * outputs_{organ}_{first_author}_{year}/
 *   {cluster_header}_cluster_statistics.csv
 *       - Columns:
 *           * cluster: Cluster name (drives parallelization)
 *           * n_cells: Number of cells in cluster
 *           * percentage: Percentage of total cells
 *       - Rows: One per cluster (sorted by cluster name)
 *
 * Example Output:
 * ---------------
 * cluster,n_cells,percentage
 * Podocyte,1234,5.2
 * Proximal Tubule,3456,14.6
 * Loop of Henle,2345,9.9
 * ...
 *
 * Resources:
 * ----------
 * - CPUs: 2 (light computation)
 * - Memory: 4-8 GB
 * - Time: ~1-2 minutes
 *
 * Notes:
 * ------
 * This is the second step in the NSForest workflow. The statistics are used
 * both for quality control and as the basis for parallelizing prep_medians.
 * The h5ad file continues through the pipeline for subsequent analyses.
 */
process cluster_stats_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(h5ad)
    
    output:
    tuple val(meta), 
          path(h5ad),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_cluster_statistics.csv"),
          emit: stats
    
    script:
    """
    nsforest-cli cluster-stats \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
