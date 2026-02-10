#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { dendrogram_process }          from '../modules/nsforest/dendrogram.nf'
include { cluster_stats_process }       from '../modules/nsforest/cluster_stats.nf'
include { prep_medians_process }        from '../modules/nsforest/prep_medians.nf'
include { merge_medians_process }       from '../modules/nsforest/merge_medians.nf'
include { prep_binary_scores_process }  from '../modules/nsforest/prep_binary_scores.nf'

// Default to test data location
params.datasets_csv = "${projectDir}/../cell-kn/data/test/kidney/cellxgene-harvester/cellxgene-harvester.csv"
params.organ = 'kidney'
params.tissue = 'kidney'
params.outdir = './test_results'
params.publish_mode = 'copy'

workflow {

  // Validate required parameters
  if (!params.datasets_csv) {
      log.error "ERROR: --datasets_csv is required"
      exit 1
  }

  if (!params.organ) {
      log.error "ERROR: --organ is required"
      exit 1
  }

  if (!params.tissue) {
      log.error "ERROR: --tissue is required"
      exit 1
  }
  
  // Read CSV and create meta map (nf-core best practice)
  csv_rows_ch = Channel
      .fromPath(params.datasets_csv)
      .ifEmpty { exit 1, "Cannot find required datasets input file : ${params.datasets_csv}" }
      .splitCsv(header: true, sep: ',')
      .map { row ->
          // Get CSV file directory to resolve relative paths
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

  // Step 1: Dendrogram
  dendrogram_output_ch = dendrogram_process(csv_rows_ch)
  
  // Step 2: Cluster stats
  cluster_stats_output_ch = cluster_stats_process(csv_rows_ch)
  
  // Step 3: SCATTER - Read clusters from stats CSV
  scattered_clusters_ch = cluster_stats_output_ch
      .flatMap { meta, h5ad, stats_csv ->
          stats_csv.splitCsv(header: true).collect { row ->
              tuple(meta, h5ad, row.cluster)
          }
      }
  
  // PROCESS: Prep medians in parallel (one per cluster)
  prep_medians_output_ch = prep_medians_process(scattered_clusters_ch)
  
  // Step 4: GATHER: Group by meta (dataset)
  gathered_medians_ch = prep_medians_output_ch.groupTuple()
  merged_medians_ch = merge_medians_process(gathered_medians_ch)
  
  // Step 5: Binary scores
  binary_scores_ch = prep_binary_scores_process(merged_medians_ch)
  
  // Debug output
  binary_scores_ch.view { meta, scores_csv ->
      "Binary scores: ${scores_csv}"
  }
}
