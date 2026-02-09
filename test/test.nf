#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { dendrogram_process }    from '../modules/nsforest/dendrogram.nf'
include { cluster_stats_process } from '../modules/nsforest/cluster_stats.nf'
include { prep_medians_process }  from '../modules/nsforest/prep_medians.nf'

params.publish_mode = 'copy'

workflow {

  // the input file, datasets_csv, for this workflow comes from the output of cellxgene-harvester https://github.com/NIH-NLM/cellxgene-harvester
  // the input parameter, organ, comes from user and will reflect be used in the filenames, it tends to be specified atm by the user, consistent with harvester
  // the input parameter, tissue, also comes from the user and should match the tissue names that are the same used by the cellxgene-harvester

  // validate required parameters
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
  
  def csv_rows_ch =
      Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find required datasets input file : ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
	    def author_cell_type_ch     = row.author_cell_type
	    def embedding_ch            = row.embedding
	    def first_author_ch         = row.first_author
	    def year_ch                 = row.year
            def h5ad_ch                 = file(row.h5ad_file)
            def disease_ch              = row.disease
            def filter_ch               = row.filter_normal
	    
        // final array for the channel
        [ author_cell_type_ch,
	  embedding_ch,
	  first_author_ch,
	  year_ch,
          h5ad_ch,
	  disease_ch,
	  filter_ch,
        ]

        organ_ch = params.organ

        tissue_ch = params.tissue
  
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
    
    // GATHER: Group by meta (dataset)
    gathered_medians_ch = prep_medians_output_ch.groupTuple()
}
