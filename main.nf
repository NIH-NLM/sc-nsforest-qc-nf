#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cluster_stats_process }          from './modules/nsforest/cluster_stats.nf'
include { compute_silhouette_process }     from './modules/scsilhouette/compute_silhouette.nf'
include { dendrogram_process }             from './modules/nsforest/dendrogram.nf'
include { download_h5ad_process }          from './modules/nsforest/download_h5ad.nf'
include { filter_adata_process }           from './modules/nsforest/filter_adata.nf'
include { merge_nsforest_results_process } from './modules/nsforest/merge_nsforest_results.nf'
include { prep_medians_process }           from './modules/nsforest/prep_medians.nf'
include { prep_binary_scores_process }     from './modules/nsforest/prep_binary_scores.nf'
include { plot_histograms_process }        from './modules/nsforest/plot_histograms.nf'
include { plots_process }                  from './modules/nsforest/plots.nf'
include { publish_results_process }        from './modules/publish/publish_results.nf'
include { run_nsforest_process }           from './modules/nsforest/run_nsforest.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_dotplot_process }            from './modules/scsilhouette/viz_dotplot.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'

params.batch_size       = 10
params.datasets_csv     = null
params.organ            = null
params.uberon_json      = null
params.disease_json     = null
params.hsapdv_json      = null
params.min_cluster_size = 5
params.outdir           = './results'
params.publish_mode     = 'copy'

workflow {
    log.info "workflow.workDir    : ${workflow.workDir}"
    log.info "workflow.launchDir  : ${workflow.launchDir}"
    log.info "workflow.projectDir : ${workflow.projectDir}"
    log.info "workflow.runName    : ${workflow.runName}"
    log.info "params.outdir       : ${params.outdir}"

    if (!params.datasets_csv) { log.error "ERROR: --datasets_csv is required"; exit 1 }
    if (!params.organ)        { log.error "ERROR: --organ is required";         exit 1 }
    if (!params.uberon_json)  { log.error "ERROR: --uberon_json is required";   exit 1 }
    if (!params.disease_json) { log.error "ERROR: --disease_json is required";  exit 1 }
    if (!params.hsapdv_json)  { log.error "ERROR: --hsapdv_json is required";   exit 1 }

    uberon_ch  = Channel.value(file(params.uberon_json))
    disease_ch = Channel.value(file(params.disease_json))
    hsapdv_ch  = Channel.value(file(params.hsapdv_json))

    csv_rows_ch = Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find datasets CSV: ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .filter { row ->
            def ref = row.reference?.trim()?.toLowerCase()
            if (ref in ['exclude', 'delete', 'merge', 'question']) {
                log.info "Skipping ${row.first_author} ${row.year} — reference='${row.reference}'"
                return false
            }
            if (!(ref in ['yes', 'no', 'unk'])) {
                log.warn "Skipping ${row.first_author} ${row.year} — unrecognised reference value '${row.reference}'"
                return false
            }
            return true
        }
        .map { row ->
            def meta = [
                organ:                              params.organ,
                first_author:                       row.first_author,
                year:                               row.year,
                author_cell_type:                   row.author_cell_type,
                embedding:                          row.embedding,
                disease:                            row.disease,
                filter:                             row.filter_normal,
                doi:                                row.doi,
                collection_name:                    row.collection_name,
                dataset_title:                      row.dataset_title,
                journal:                            row.journal,
                collection_url:                     row.collection_url,
                explorer_url:                       row.explorer_url,
                h5ad_url:                           row.h5ad_url,
                tissue_ontology_term_id:            row.tissue_ontology_term_id,
                disease_ontology_term_id:           row.disease_ontology_term_id,
                development_stage_ontology_term_id: row.development_stage_ontology_term_id,
            ]
            tuple(meta, row.h5ad_url)
        }

    // Step 0a: Download h5ad from CellxGene URL
    downloaded_ch = download_h5ad_process(csv_rows_ch)

    // Step 0b: Filter — tissue + disease + age using per-row ontology term IDs
    filter_output_ch = filter_adata_process(
        downloaded_ch.h5ad
            .combine(uberon_ch)
            .combine(disease_ch)
            .combine(hsapdv_ch)
    )

    // Convenience: filtered h5ad only channel
    filtered_h5ad_ch = filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }

    // Step 1: Dendrogram
    dendrogram_output_ch = dendrogram_process(filtered_h5ad_ch)

    // Step 1b: Cluster statistics
    cluster_stats_output_ch = cluster_stats_process(filtered_h5ad_ch)
    
    // Step 2a: Prep medians — runs once per dataset on full filtered h5ad
    prep_medians_output_ch = prep_medians_process(filtered_h5ad_ch)

    // Step 2b: Prep binary scores — runs once per dataset on full filtered h5ad
    prep_binary_scores_output_ch = prep_binary_scores_process(filtered_h5ad_ch)

    // Step 3: Plot histograms
    plot_histograms_process(
        prep_medians_output_ch.complete
            .map { meta, medians_csv, medians_pkl -> tuple(meta, medians_csv) }
            .join(
                prep_binary_scores_output_ch.complete
                    .map { meta, binary_csv, binary_pkl -> tuple(meta, binary_csv) }
            )
            .map { meta, medians_csv, binary_csv -> tuple(meta, medians_csv, binary_csv) }
    )

    // Step 4: Scatter run_nsforest by cluster batch
    def batchSize = params.batch_size ?: 10
    nsforest_input_ch = filtered_h5ad_ch
        .join(prep_medians_output_ch.complete.map { meta, medians_csv, medians_pkl -> tuple(meta, medians_csv) })
        .join(prep_binary_scores_output_ch.complete.map { meta, binary_csv, binary_pkl -> tuple(meta, binary_csv) })
        .join(dendrogram_output_ch.stats.map { meta, h5ad, cluster_order_csv -> tuple(meta, cluster_order_csv) })
        .flatMap { meta, h5ad, medians_csv, binary_csv, cluster_order_csv ->
            def clusters = cluster_order_csv
                .splitCsv(header: true)
                .collect { it.cluster_order }
            clusters.collate(batchSize).collect { batch ->
                tuple(meta, h5ad, medians_csv, binary_csv, batch.join(','))
            }
        }

    nsforest_output_ch = run_nsforest_process(nsforest_input_ch)

    // Step 5: Merge NSForest results
    merged_nsforest_ch = merge_nsforest_results_process(
        nsforest_output_ch.partial.groupTuple()
    )

    // Step 6: Plots
    plots_process(
        filtered_h5ad_ch
            .join(
                merged_nsforest_ch.complete.map { meta, results_csv, results_pkl -> tuple(meta, results_csv) }
            )
            .map { meta, h5ad, results_csv -> tuple(meta, h5ad, results_csv) }
    )

    // Step 7: Compute silhouette
    silhouette_output_ch = compute_silhouette_process(filtered_h5ad_ch)

    // Step 8a: viz_summary
    viz_summary_process(
        silhouette_output_ch.results
            .map { meta, files ->
                def flist      = files instanceof List ? files : [files]
                def scores     = flist.find { it.name.endsWith('_silhouette_scores.csv') }
                def summary    = flist.find { it.name.endsWith('_cluster_summary.csv') }
                def annotation = flist.find { it.name.endsWith('_annotation.json') }
                tuple(meta, scores, summary, annotation)
            }
            .join(
                merged_nsforest_ch.complete.map { meta, results_csv, results_pkl -> tuple(meta, results_csv) },
                remainder: true
            )
            .map { meta, scores, summary, annotation, nsforest_csv ->
                tuple(meta, scores, summary, annotation, nsforest_csv ?: file('NO_FILE'))
            }
    )

    // Step 8b: viz_distribution
    viz_distribution_process(
        silhouette_output_ch.results.map { meta, files ->
            def flist      = files instanceof List ? files : [files]
            def scores     = flist.find { it.name.endsWith('_silhouette_scores.csv') }
            def summary    = flist.find { it.name.endsWith('_cluster_summary.csv') }
            def annotation = flist.find { it.name.endsWith('_annotation.json') }
            tuple(meta, scores, summary, annotation)
        }
    )

    // Step 8c: viz_dotplot
    viz_dotplot_process(filtered_h5ad_ch)

    // Step 9: Publish
    if (params.github_token) {
        all_files_ch = Channel
            .empty()
            .mix(
	        dendrogram_output_ch.stats.map              { meta, files -> tuple(meta, files) },
	        dendrogram_output_ch.results.map            { meta, files -> tuple(meta, files) },
		cluster_stats_process.results.map           { meta, files -> tuple(meta, files) },
		filter_adata_process.results.map            { meta, files -> tuple(meta, files) },
		merge_results_process.complete.map          { meta, files -> tuple(meta, files) },
		plot_histograms.histograms.map              { meta, files -> tuple(meta, files) },
		plots_process.plots.map                     { meta, files -> tuple(meta, files) },
		prep_medians_process.complete.map           { meta, files -> tuple(meta, files) },
		prep_binary_scores_process.complete.map     { meta, files -> tuple(meta, files) },
                merge_nsforest_results_process.complete.map { meta, files -> tuple(meta, files) },
		run_nsforest_process.partial.map            { meta, files -> tuple(meta, files) },
                compute_silhouette_process.results.map      { meta, files -> tuple(meta, files) },
                viz_dotplot_process.plots.map               { meta, files -> tuple(meta, files) },
                viz_distribution_process.plots.map          { meta, files -> tuple(meta, files) },
                viz_summary_process.plots.map               { meta, files -> tuple(meta, files) },
            )
            .map { meta, file_lists ->
                tuple(meta[0], file_lists.flatten())
            }
        publish_results_process(all_files_ch)
    } else {
        log.warn "WARNING: --github_token not set — skipping publish step"
    }
}
