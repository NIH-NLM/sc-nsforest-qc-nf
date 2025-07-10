#!/usr/bin/env nextflow

include { nsforest_process }  from './modules/nsforest_process.nf'

workflow {

  def csv_rows_ch =
      Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find required datasets input file : ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def h5ad_ch             = file(row.h5ad_file)
            def label_key_ch        = row.label_key
            def embedding_key_ch    = row.embedding_key
            def organism_ch         = row.organism
            def disease_ch          = row.disease
            def tissue_ch           = row.tissue
	    def author_ch           = row.author
	    def publication_date_ch = row.publication_date
	    def publication_ch      = row.publication
            def cell_count_ch       = row.cell_count
	    
        // final array for the channel
        [ h5ad_ch, label_key_ch, embedding_key_ch , organism_ch, disease_ch, tissue_ch,
	  author_ch, publication_date_ch, publication_ch, cell_count_ch ]
      }
	
   nsforest_process (
        csv_rows_ch )

   nsforest_process.out.nsforest_process_output_ch
}
