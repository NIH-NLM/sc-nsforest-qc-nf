#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * NSForest Workflow - Test Pipeline
 *
 * Complete workflow from filtering through plotting with parallelization.
 * Uses tuple passing to carry adata through pipeline stages.
 */

// Include all modules
include { filter_adata_process }         from '../modules/nsforest/filter_adata.nf'
include { prep_medians_process }         from '../modules/nsforest/prep_medians.nf'
include { merge_medians_process }        from '../modules/nsforest/merge_medians.nf'
include { prep_binary_scores_process }   from '../modules/nsforest/prep_binary_scores.nf'
include { merge_binary_scores_process }  from '../modules/nsforest/merge_binary_scores.nf'
include { plot_histograms_process }      from '../modules/nsforest/plot_histograms.nf'
include { run_nsforest_process }         from '../modules/nsforest/run_nsforest.nf'
include { merge_nsforest_results_process } from '../modules/nsforest/merge_nsforest_results.nf'
include { plots_process }                from '../modules/nsforest/plots.nf'

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
    
    // Step 1: Filter adata
    filter_output_ch = filter_adata_process(csv_rows_ch)
    
    filter_output_ch.results.view { meta, adata_filtered, stats ->
        "Step 1 - Filtered: ${adata_filtered}"
    }
    
    // Extract filtered h5ad and cluster sizes
    filtered_data_ch = filter_output_ch.results
        .map { meta, adata_filtered, stats_files ->
            def cluster_sizes_file = stats_files.find { 
                it.name.contains('cluster_sizes.csv') && !it.name.contains('before_filter') 
            }
            tuple(meta, adata_filtered, cluster_sizes_file)
        }
    
    // Step 2: SCATTER prep_medians by cluster
    scattered_medians_ch = filtered_data_ch
        .flatMap { meta, adata_filtered, cluster_sizes_csv ->
            cluster_sizes_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_filtered, row.cluster)
            }
        }
    
    scattered_medians_ch.view { meta, adata, cluster ->
        "Step 2 - Scatter medians for cluster: ${cluster}"
    }
    
    prep_medians_output_ch = prep_medians_process(scattered_medians_ch)
    
    // Step 3: GATHER merge_medians
    gathered_medians_ch = prep_medians_output_ch.partial
        .groupTuple()
    
    gathered_medians_ch.view { meta, adata_preps, csvs ->
        "Step 3 - Gathered ${csvs.size()} median files"
    }
    
    merged_medians_ch = merge_medians_process(gathered_medians_ch)
    
    merged_medians_ch.complete.view { meta, adata_prep, files ->
        "Step 3 - Merged medians: ${adata_prep}, files: ${files}"
    }
    
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
    
    scattered_binary_ch.view { meta, adata_prep, cluster ->
        "Step 4 - Scatter binary scores for cluster: ${cluster}"
    }
    
    prep_binary_output_ch = prep_binary_scores_process(scattered_binary_ch)
    
    // Step 5: GATHER merge_binary_scores
    gathered_binary_ch = prep_binary_output_ch.partial
        .groupTuple()
    
    gathered_binary_ch.view { meta, adata_preps, csvs ->
        "Step 5 - Gathered ${csvs.size()} binary score files"
    }
    
    merged_binary_ch = merge_binary_scores_process(gathered_binary_ch)
    
    merged_binary_ch.complete.view { meta, adata_prep, files ->
        "Step 5 - Merged binary scores: ${adata_prep}, files: ${files}"
    }
    
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
    
    histogram_input_ch.view { meta, medians_csv, binary_csv ->
        "Step 6 - Plotting histograms: ${medians_csv}, ${binary_csv}"
    }
    
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
    
    scattered_nsforest_ch.view { meta, adata_prep, cluster ->
        "Step 7 - Scatter NSForest for cluster: ${cluster}"
    }
    
    nsforest_output_ch = run_nsforest_process(scattered_nsforest_ch)
    
    // Step 8: GATHER merge_nsforest_results
    gathered_nsforest_ch = nsforest_output_ch.partial
        .groupTuple()
    
    gathered_nsforest_ch.view { meta, csvs ->
        "Step 8 - Gathered ${csvs.size()} NSForest result files"
    }
    
    merged_nsforest_ch = merge_nsforest_results_process(gathered_nsforest_ch)
    
    merged_nsforest_ch.complete.view { meta, files ->
        "Step 8 - NSForest results merged: ${files}"
    }
    
    // Step 9: Plots (needs original filtered h5ad + results)
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
    
    plots_input_ch.view { meta, adata, results_csv ->
        "Step 9 - Creating plots with: ${adata}, ${results_csv}"
    }
    
    plots_process(plots_input_ch)
}
