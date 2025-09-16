#!/usr/bin/env nextflow

include { build_symbol_map_process }    from './modules/build_symbol_map.nf'
include { run_dendrogramplot_process }  from './modules/dendrogram.nf'
include { run_dotplot_process }         from './modules/dotplot.nf'
include { filter_condition_process }    from './modules/filter_condition.nf'
include { filter_tissue_process }       from './modules/filter_tissue.nf'
include { run_matrixplot_process }      from './modules/matrixplot.nf'
include { prep_binary_scores_process }  from './modules/prep_binary_scores.nf' 
include { prep_medians_process }        from './modules/prep_medians.nf'
include { run_nsforest_process }        from './modules/run_nsforest.nf'
include { sanitize_labels_process }     from './modules/sanitize_labels.nf'
include { symbolize_genes_process }     from './modules/symbolize_genes.nf'
include { run_violinplot_process }      from './modules/violinplot.nf'

workflow {

  def csv_rows_ch =
      Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find required datasets input file : ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def h5ad_ch                 = file(row.h5ad_file)
            def base_ch                 = h5ad_ch.getSimpleName().replaceFirst(/\.h5ad$/, '')
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
	  tissue_ch, author_ch, publication_date_ch, publication_ch, cell_count_ch, base_ch ]
      }

      sanitize_output_ch       = sanitize_labels_process (
         csv_rows_ch )

      filter_disease_output_ch = filter_condition_process (
         sanitize_output_ch )

      filter_tissue_output_ch  = filter_tissue_process (
         filter_disease_output_ch )

      prep_medians_output_ch   = prep_medians_process (
         filter_tissue_output_ch )

      binary_scores_output_ch = prep_binary_scores_process (
         prep_medians_output_ch )

      nsforest_output_ch = run_nsforest_process (
         binary_scores_output_ch )

      symbol_map_csv_ch = build_symbol_map_process (
         nsforest_output_ch )

      symbolize_genes_ch = symbolize_genes_process (
          symbol_map_csv_ch )

      run_dendrogram_output_ch = run_dendrogramplot_process (
         symbolize_genes_ch )

      run_dotplot_process (
         run_dendrogram_output_ch )

      run_matrixplot_process (
         run_dendrogram_output_ch )

     run_violinplot_process (
         run_dendrogram_output_ch )

}
