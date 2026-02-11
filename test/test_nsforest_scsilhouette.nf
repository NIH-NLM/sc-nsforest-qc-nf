#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * NSForest + scsilhouette QC Pipeline - Complete Workflow
 *
 * Data flow:
 * 1. Filter adata by disease/tissue/min_cluster_size → adata_filtered.h5ad
 * 2. NSForest pipeline uses adata_filtered → marker discovery
 * 3. scsilhouette uses adata_filtered → clustering QC
 * 4. Integrate NSForest + scsilhouette results → combined QC visualization
 * 5. Generate dataset summary statistics → median-of-medians
 */

// Include NSForest modules
include { filter_adata_process }         from '../modules/nsforest/filter_adata.nf'
include { prep_medians_process }         from '../modules/nsforest/prep_medians.nf'
include { merge_medians_process }        from '../modules/nsforest/merge_medians.nf'
include { prep_binary_scores_process }   from '../modules/nsforest/prep_binary_scores.nf'
include { merge_binary_scores_process }  from '../modules/nsforest/merge_binary_scores.nf'
include { plot_histograms_process }      from '../modules/nsforest/plot_histograms.nf'
include { run_nsforest_process }         from '../modules/nsforest/run_nsforest.nf'
include { merge_nsforest_results_process } from '../modules/nsforest/merge_nsforest_results.nf'
include { plots_process }                from '../modules/nsforest/plots.nf'

// Include scsilhouette modules
include { compute_silhouette_process }   from '../modules/scsilhouette/compute_silhouette.nf'
include { viz_summary_process }          from '../modules/scsilhouette/viz_summary.nf'
include { viz_dotplot_process }          from '../modules/scsilhouette/viz_dotplot.nf'
include { viz_distribution_process }     from '../modules/scsilhouette/viz_distribution.nf'
include { compute_summary_stats_process } from '../modules/scsilhouette/compute_summary_stats.nf'

// Parameters
params.datasets_csv = null
params.organ = null
params.tissue = null
params.outdir = './results'
params.publish_mode = 'copy'

workflow {
    // Validate required parameters
    if (!params.datasets_csv) {
        error "Missing required parameter: --datasets_csv"
    }
    if (!params.organ) {
        error "Missing required parameter: --organ"
    }
    if (!params.tissue) {
        error "Missing required parameter: --tissue"
    }
    
    // Read CSV and create input channel with meta map
    csv_rows_ch = Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find required datasets input file: ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def csv_dir = file(params.datasets_csv).parent
            def h5ad_path = row.h5ad_url.startsWith('/') ? 
                            file(row.h5ad_url) : 
                            file("${csv_dir}/${row.h5ad_url}")
            
            def meta = [
                organ           : params.organ,
                tissue          : params.tissue,
                author_cell_type: row.author_cell_type,
                embedding       : row.embedding,
                first_author    : row.first_author,
                year            : row.year,
                disease         : row.disease,
                filter          : row.filter_normal
            ]
            tuple(meta, h5ad_path)
        }
    
    // ========================================================================
    // SHARED FILTERING STEP
    // Both NSForest and scsilhouette use the same filtered adata
    // ========================================================================
    
    // Step 1: Filter adata by disease, tissue, min cluster size
    // Output: adata_filtered.h5ad (used by BOTH NSForest and scsilhouette)
    filter_output_ch = filter_adata_process(csv_rows_ch)
    
    filter_output_ch.results.view { meta, adata_filtered, stats ->
        "SHARED Step 1 - Filtered adata (used by both NSForest & scsilhouette): ${adata_filtered}"
    }
    
    filtered_data_ch = filter_output_ch.results
        .map { meta, adata_filtered, stats_files ->
            def cluster_sizes_file = stats_files.find { 
                it.name.contains('cluster_sizes.csv') && !it.name.contains('before_filter') 
            }
            tuple(meta, adata_filtered, cluster_sizes_file)
        }
    
    // ========================================================================
    // NSFOREST WORKFLOW - Uses adata_filtered.h5ad
    // ========================================================================
    
    // Step 2: SCATTER prep_medians by cluster
    scattered_medians_ch = filtered_data_ch
        .flatMap { meta, adata_filtered, cluster_sizes_csv ->
            cluster_sizes_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_filtered, row.cluster)
            }
        }
    
    prep_medians_output_ch = prep_medians_process(scattered_medians_ch)
    
    // Step 3: GATHER merge_medians
    gathered_medians_ch = prep_medians_output_ch.partial.groupTuple()
    merged_medians_ch = merge_medians_process(gathered_medians_ch)
    
    // Step 4: SCATTER prep_binary_scores by cluster
    scattered_binary_ch = merged_medians_ch.complete
        .map { meta, adata_prep, medians_files ->
            tuple(meta, adata_prep)
        }
        .combine(
            filtered_data_ch.map { meta, adata_filtered, cluster_sizes_csv ->
                tuple(meta, cluster_sizes_csv)
            },
            by: 0
        )
        .flatMap { meta, adata_prep, cluster_sizes_csv ->
            cluster_sizes_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }
    
    prep_binary_output_ch = prep_binary_scores_process(scattered_binary_ch)
    
    // Step 5: GATHER merge_binary_scores
    gathered_binary_ch = prep_binary_output_ch.partial.groupTuple()
    merged_binary_ch = merge_binary_scores_process(gathered_binary_ch)
    
    // Step 6: Plot histograms
    histogram_input_ch = merged_medians_ch.complete
        .map { meta, adata_prep, medians_files ->
            def medians_csv = medians_files.find { it.name.endsWith('.csv') }
            tuple(meta, medians_csv)
        }
        .join(
            merged_binary_ch.complete.map { meta, adata_prep, binary_files ->
                def binary_csv = binary_files.find { it.name.endsWith('.csv') }
                tuple(meta, binary_csv)
            }
        )
    
    plot_histograms_process(histogram_input_ch)
    
    // Step 7: SCATTER run_nsforest by cluster
    scattered_nsforest_ch = merged_binary_ch.complete
        .map { meta, adata_prep, binary_files ->
            tuple(meta, adata_prep)
        }
        .combine(
            filtered_data_ch.map { meta, adata_filtered, cluster_sizes_csv ->
                tuple(meta, cluster_sizes_csv)
            },
            by: 0
        )
        .flatMap { meta, adata_prep, cluster_sizes_csv ->
            cluster_sizes_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }
    
    nsforest_output_ch = run_nsforest_process(scattered_nsforest_ch)
    
    // Step 8: GATHER merge_nsforest_results
    gathered_nsforest_ch = nsforest_output_ch.partial.groupTuple()
    merged_nsforest_ch = merge_nsforest_results_process(gathered_nsforest_ch)
    
    merged_nsforest_ch.complete.view { meta, files ->
        "NSForest Step 8 - Results merged: ${files}"
    }
    
    // Step 9: NSForest plots (uses original adata_filtered for expression plots)
    plots_input_ch = filtered_data_ch
        .map { meta, adata_filtered, cluster_sizes_csv ->
            tuple(meta, adata_filtered)
        }
        .join(
            merged_nsforest_ch.complete.map { meta, result_files ->
                def results_csv = result_files.find { it.name.endsWith('.csv') }
                tuple(meta, results_csv)
            }
        )
    
    plots_process(plots_input_ch)
    
    // ========================================================================
    // SCSILHOUETTE WORKFLOW - Uses adata_filtered.h5ad (same filtered data)
    // ========================================================================
    
    // Step 10: Compute silhouette scores on FILTERED adata
    silhouette_input_ch = filtered_data_ch
        .map { meta, adata_filtered, cluster_sizes_csv ->
            tuple(meta, adata_filtered)
        }
    
    silhouette_output_ch = compute_silhouette_process(silhouette_input_ch)
    
    silhouette_output_ch.results.view { meta, scores, summary, annotation ->
        "scsilhouette Step 10 - Computed scores on FILTERED adata: ${scores}"
    }
    
    // Step 11: Silhouette embedding dotplot (uses FILTERED adata)
    dotplot_input_ch = filtered_data_ch
        .map { meta, adata_filtered, cluster_sizes_csv ->
            tuple(meta, adata_filtered)
        }
    
    viz_dotplot_process(dotplot_input_ch)
    
    // Step 12: Distribution plots
    distribution_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
    
    viz_distribution_process(distribution_input_ch)
    
    // Step 13: Summary visualization WITH NSForest F-scores integration
    summary_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
        .join(
            merged_nsforest_ch.complete.map { meta, result_files ->
                def results_csv = result_files.find { it.name.endsWith('.csv') }
                tuple(meta, results_csv)
            }
        )
    
    summary_input_ch.view { meta, scores, summary, annotation, nsforest_results ->
        "scsilhouette Step 13 - Creating integrated summary with NSForest results: ${nsforest_results}"
    }
    
    viz_summary_process(summary_input_ch)
    
    // ========================================================================
    // FINAL SUMMARY STATISTICS
    // Dataset-level summary: median-of-medians, cluster quality distribution
    // ========================================================================
    
    // Step 14: Compute dataset summary statistics
    dataset_summary_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
        .join(
            merged_nsforest_ch.complete.map { meta, result_files ->
                def results_csv = result_files.find { it.name.endsWith('.csv') }
                tuple(meta, results_csv)
            }
        )
    
    compute_summary_stats_process(dataset_summary_input_ch)
    
    dataset_summary_input_ch.view { meta, scores, summary, annotation, nsforest ->
        "FINAL Step 14 - Computing dataset summary (median-of-medians)"
    }
}
