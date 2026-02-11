#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * scsilhouette Standalone Test - Using NSForest Outputs
 *
 * Tests scsilhouette workflow using existing filtered adata and NSForest results.
 * No filtering needed - uses adata_filtered.h5ad from previous NSForest run.
 */

// Include scsilhouette modules
include { compute_silhouette_process }   from '../modules/scsilhouette/compute_silhouette.nf'
include { viz_summary_process }          from '../modules/scsilhouette/viz_summary.nf'
include { viz_dotplot_process }          from '../modules/scsilhouette/viz_dotplot.nf'
include { viz_distribution_process }     from '../modules/scsilhouette/viz_distribution.nf'
include { compute_summary_stats_process } from '../modules/scsilhouette/compute_summary_stats.nf'

// Parameters
params.nsforest_results_dir = null  // Path to outputs_kidney_Lake_2023/
params.organ = null
params.first_author = null
params.year = null
params.author_cell_type = null
params.embedding = null
params.outdir = './results'
params.publish_mode = 'copy'

workflow {
    // Validate required parameters
    if (!params.nsforest_results_dir) {
        error "Missing required parameter: --nsforest_results_dir (e.g., results/outputs_kidney_Lake_2023)"
    }
    if (!params.organ) {
        error "Missing required parameter: --organ"
    }
    if (!params.first_author) {
        error "Missing required parameter: --first_author"
    }
    if (!params.year) {
        error "Missing required parameter: --year"
    }
    if (!params.author_cell_type) {
        error "Missing required parameter: --author_cell_type"
    }
    if (!params.embedding) {
        error "Missing required parameter: --embedding (e.g., X_umap)"
    }
    
    // Create meta map
    def meta = [
        organ           : params.organ,
        tissue          : params.organ,  // Assume same
        author_cell_type: params.author_cell_type,
        embedding       : params.embedding,
        first_author    : params.first_author,
        year            : params.year,
        disease         : null,
        filter          : null
    ]
    
    // Locate input files from NSForest results
    def results_dir = file(params.nsforest_results_dir)
    def adata_filtered = file("${results_dir}/adata_filtered.h5ad")
    def nsforest_results = file("${results_dir}/${params.author_cell_type}_results.csv")
    
    // Verify files exist
    if (!adata_filtered.exists()) {
        error "Cannot find adata_filtered.h5ad in ${params.nsforest_results_dir}"
    }
    if (!nsforest_results.exists()) {
        error "Cannot find ${params.author_cell_type}_results.csv in ${params.nsforest_results_dir}"
    }
    
    log.info """
    ========================================
    scsilhouette Standalone Test
    ========================================
    Results directory: ${params.nsforest_results_dir}
    Filtered adata   : ${adata_filtered}
    NSForest results : ${nsforest_results}
    Embedding        : ${params.embedding}
    Output directory : ${params.outdir}
    ========================================
    """.stripIndent()
    
    // Create input channels
    adata_ch = Channel.value(tuple(meta, adata_filtered))
    
    // ========================================================================
    // SCSILHOUETTE WORKFLOW
    // ========================================================================
    
    // Step 1: Compute silhouette scores
    silhouette_output_ch = compute_silhouette_process(adata_ch)
    
    silhouette_output_ch.results.view { meta, scores, summary, annotation ->
        "Step 1 - Computed silhouette scores: ${scores}"
    }
    
    // Step 2: Embedding dotplot
    viz_dotplot_process(adata_ch)
    
    // Step 3: Distribution plots
    distribution_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
    
    viz_distribution_process(distribution_input_ch)
    
    // Step 4: Summary visualization with NSForest F-scores
    nsforest_ch = Channel.value(tuple(meta, nsforest_results))
    
    summary_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
        .combine(nsforest_ch, by: 0)
    
    summary_input_ch.view { meta, scores, summary, annotation, nsforest ->
        "Step 4 - Creating integrated summary with NSForest results"
    }
    
    viz_summary_process(summary_input_ch)
    
    // Step 5: Dataset summary statistics
    compute_summary_stats_process(summary_input_ch)
}
