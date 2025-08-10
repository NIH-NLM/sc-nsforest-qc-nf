#!/usr/bin/env nextflow

include { prep_medians_process }       from './modules/prep_medians.nf'
include { prep_binary_scores_process } from './modules/prep_binary_scores.nf' 
include { run_nsforest_process }       from './modules/run_nsforest.nf'

workflow {

  def csv_rows_ch =
      Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find required datasets input file : ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def h5ad_ch                 = file(row.h5ad_file)
            def label_key_ch            = row.label_key
            def embedding_key_ch        = row.embedding_key
            def organism_ch             = row.organism
            def disease_ch              = row.disease
            def filter_ch               = row.filter_normal
            def metric_ch               = row.metric
            def save_scores_ch          = row.save_scores
            def save_cluster_summary_ch = row.save_cluster_summary
            def save_annotation_ch      = row.save_annotation
            def tissue_ch               = row.tissue
            def author_ch               = row.author
            def publication_date_ch     = row.publication_date
	    def publication_ch          = row.publication
            def cell_count_ch           = row.cell_count
	    
        // final array for the channel
        [ h5ad_ch, label_key_ch, embedding_key_ch , organism_ch, disease_ch,
	  filter_ch, metric_ch, save_scores_ch, save_cluster_summary_ch, save_annotation_ch,
	  tissue_ch, author_ch, publication_date_ch, publication_ch, cell_count_ch ]
      }
	
      prep_medians_process (
        csv_rows_ch )

      prep_binary_scores_process (
        prep_medians_process.out.prep_medians_output_ch )

      run_nsforest_process (
        prep_binary_scores_process.out.prep_binary_scores_output_ch )


}
