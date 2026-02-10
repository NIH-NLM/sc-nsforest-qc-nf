#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// NSForest processes
include { dendrogram_process }                     from './modules/nsforest/dendrogram.nf'
include { cluster_stats_process }                  from './modules/nsforest/cluster_stats.nf'
include { prep_medians_process }                   from './modules/nsforest/prep_medians.nf'
include { prep_binary_scores_process }             from './modules/nsforest/prep_binary_scores.nf'
include { plot_histograms_process }                from './modules/nsforest/plot_histograms.nf'
include { run_nsforest_process }                   from './modules/nsforest/run_nsforest.nf'
include { plot_boxplots_process }                  from './modules/nsforest/plot_boxplots.nf'
include { plot_scatter_process }                   from './modules/nsforest/plot_scatter.nf'
include { gene_mapping_process }                   from './modules/nsforest/gene_mapping.nf'
include { plot_expression_dotplot_process }        from './modules/nsforest/plot_expression_dotplot.nf'
include { plot_expression_stacked_violin_process } from './modules/nsforest/plot_matrixplot.nf'

// scsilhouette processes
include { compute_silhouette_process }             from './modules/scsilhouette/compute_scsilhouette.nf'
include { plot_qc_summary_process }                from './modules/scsilhouette/plot_qc_summary.nf'
include { plot_dotplot_process }                   from './modules/scsilhouette/plot_dotplot.nf'
include { plot_distribution_process }              from './modules/scsilhouette/plot_distribution.nf'

params.datasets_csv = null
params.organ = null
params.tissue = null
params.outdir = './results'
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
          def meta = [
              organ: params.organ,
              first_author: row.first_author,
              year: row.year,
              author_cell_type: row.author_cell_type,
              embedding: row.embedding,
              disease: row.disease,
              filter: row.filter_normal,
              tissue: params.tissue
          ]
          tuple(meta, file(row.h5ad_file))
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
  
  // Step 4: PROCESS - Prep medians in parallel (one per cluster)
  prep_medians_output_ch = prep_medians_process(scattered_clusters_ch)
  
  // Step 5: GATHER - Group by meta (dataset)
  gathered_medians_ch = prep_medians_output_ch.groupTuple()
  
  // Continue with remaining steps...
  // prep_binary_scores_process(gathered_medians_ch)
  // etc.
}
