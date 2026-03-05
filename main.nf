#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { download_h5ad_process }                  from './modules/nsforest/download_h5ad.nf'
include { filter_adata_process }                   from './modules/nsforest/filter_adata.nf'
include { dendrogram_process }                     from './modules/nsforest/dendrogram.nf'
include { prep_medians_process }                   from './modules/nsforest/prep_medians.nf'
include { prep_binary_scores_process }             from './modules/nsforest/prep_binary_scores.nf'
include { plot_histograms_process }                from './modules/nsforest/plot_histograms.nf'
include { run_nsforest_process }                   from './modules/nsforest/run_nsforest.nf'
include { merge_nsforest_results_process }         from './modules/nsforest/merge_nsforest_results.nf'
include { plots_process }                          from './modules/nsforest/plots.nf'
include { compute_silhouette_process }             from './modules/scsilhouette/compute_silhouette.nf'
include { viz_summary_process }                    from './modules/scsilhouette/viz_summary.nf'
include { viz_dotplot_process }                    from './modules/scsilhouette/viz_dotplot.nf'
include { viz_distribution_process }               from './modules/scsilhouette/viz_distribution.nf'
include { publish_results_process }             from './modules/publish/publish_results.nf'

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
            def ref = row.reference?.trim().toLowerCase()
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
                organ:            params.organ,
                first_author:     row.first_author,
                year:             row.year,
                author_cell_type: row.author_cell_type,
                embedding:        row.embedding,
                disease:          row.disease,
                filter:           row.filter_normal,
                doi:              row.doi,
                collection_name:  row.collection_name,
                dataset_title:    row.dataset_title,
                journal:          row.journal,
                collection_url:   row.collection_url,
                explorer_url:     row.explorer_url,
                h5ad_url:         row.h5ad_url,
            ]
            tuple(meta, row.h5ad_url)
	}

    // Step 0a: Download h5ad from CellxGene URL
    downloaded_ch = download_h5ad_process(csv_rows_ch)

    // Step 0b: Filter — tissue + disease + age
    filter_output_ch = filter_adata_process(
        downloaded_ch.h5ad
            .combine(uberon_ch)
            .combine(disease_ch)
            .combine(hsapdv_ch)
    )

    // Step 1: Dendrogram — drives scatter for run_nsforest
    dendrogram_output_ch = dendrogram_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )

    // Step 2: Prep medians — runs ONCE per dataset on full filtered h5ad
    // No scatter here — global gene set must be preserved for NSForest
    prep_medians_output_ch = prep_medians_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )
    // emit: complete -> (meta, adata_prep.h5ad, medians.csv)

    // Step 3: Prep binary scores — runs ONCE per dataset
    // No scatter — needs full gene set from prep_medians
    prep_binary_scores_output_ch = prep_binary_scores_process(
        prep_medians_output_ch.complete.map { meta, adata_prep, medians -> tuple(meta, adata_prep) }
    )
    // emit: complete -> (meta, adata_prep.h5ad, binary_scores.csv)

    // Step 4: Plot histograms
    plot_histograms_process(
        prep_medians_output_ch.complete
            .map { meta, adata_prep, medians -> tuple(meta, medians) }
            .join(
                prep_binary_scores_output_ch.complete
                    .map { meta, adata_prep, binary -> tuple(meta, binary) }
            )
            .map { meta, medians_csv, binary_csv -> tuple(meta, medians_csv, binary_csv) }
    )

    // Step 5: Scatter run_nsforest by cluster batch
    // Only NSForest is scattered — safe because it reads from adata_prep independently
    def batchSize = params.batch_size ?: 10
    nsforest_input_ch = prep_binary_scores_output_ch.complete
        .map { meta, adata_prep, binary -> tuple(meta, adata_prep) }
        .combine(
            dendrogram_output_ch.stats.map { meta, h5ad, cluster_order_csv -> tuple(meta, cluster_order_csv) },
            by: 0
        )
        .flatMap { meta, adata_prep, cluster_order_csv ->
            def clusters = cluster_order_csv
                .splitCsv(header: true)
                .collect { it.cluster_order }
            clusters.collate(batchSize).collect { batch ->
                tuple(meta, adata_prep, batch.join(','))
            }
        }

    nsforest_output_ch = run_nsforest_process(nsforest_input_ch)

    // Step 6: Gather NSForest results
    merged_nsforest_ch = merge_nsforest_results_process(
        nsforest_output_ch.partial.groupTuple()
    )
    // emit: complete -> (meta, results.{csv,pkl})

    // Step 7: Plots — gene symbol mapping handled internally by nsforest-cli plots
    plots_process(
        filter_output_ch.results
            .map { meta, h5ad, stats -> tuple(meta, h5ad) }
            .join(
                merged_nsforest_ch.complete.map { meta, results_files ->
                    def results_csv = results_files instanceof List
                        ? results_files.find { it.name.endsWith('.csv') } : results_files
                    tuple(meta, results_csv)
                }
            )
            .map { meta, h5ad, results_csv -> tuple(meta, h5ad, results_csv) }
    )

    // Step 8: Compute silhouette — h5ad already filtered
    silhouette_output_ch = compute_silhouette_process(
        filter_output_ch.results
            .map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )
    
    // Step 9a: viz_summary
    viz_summary_process(
        silhouette_output_ch.results
            .map { meta, files ->
                def flist = files instanceof List ? files : [files]
                def scores     = flist.find { it.name.endsWith('_silhouette_scores.csv') }
                def summary    = flist.find { it.name.endsWith('_cluster_summary.csv') }
                def annotation = flist.find { it.name.endsWith('_annotation.json') }
                tuple(meta, scores, summary, annotation)
            }
            .join(
                merged_nsforest_ch.complete.map { meta, results_files ->
                    def results_csv = results_files instanceof List
                        ? results_files.find { it.name.endsWith('.csv') } : results_files
                    tuple(meta, results_csv)
                },
                remainder: true
            )
            .map { meta, scores, summary, annotation, nsforest_csv ->
                tuple(meta, scores, summary, annotation, nsforest_csv ?: file('NO_FILE'))
            }
    )

    // Step 9b: viz_distribution
    viz_distribution_process(
        silhouette_output_ch.results.map { meta, files ->
            def flist = files instanceof List ? files : [files]
            def scores     = flist.find { it.name.endsWith('_silhouette_scores.csv') }
            def summary    = flist.find { it.name.endsWith('_cluster_summary.csv') }
            def annotation = flist.find { it.name.endsWith('_annotation.json') }
            tuple(meta, scores, summary, annotation)
        }
    )

    // Step 9c: viz_dotplot
    viz_dotplot_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )

    // Step 10: Publish — fires once ALL datasets complete both branches
    if (params.github_token) {
        all_files_ch = plots_process.out.plots
            .map { meta, files ->
                def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                def flist = files instanceof List ? files : [files]
                flist.collect { f -> "${label}:::${f}" }
            }
            .mix(
                viz_summary_process.out.plots.map { meta, files, dataset_summary ->
                    def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                    def flist = files instanceof List ? files : [files]
                    (flist + [dataset_summary]).collect { f -> "${label}:::${f}" }
                }
            )
            .mix(
                viz_distribution_process.out.plots.map { meta, files ->
                    def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                    def flist = files instanceof List ? files : [files]
                    flist.collect { f -> "${label}:::${f}" }
                }
            )
            .mix(
                viz_dotplot_process.out.plots.map { meta, files ->
                    def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                    def flist = files instanceof List ? files : [files]
                    flist.collect { f -> "${label}:::${f}" }
                }
            )
            .mix(
                compute_silhouette_process.out.results.map { meta, scores, cluster_summary, annotation ->
                    def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                    [scores, cluster_summary, annotation].collect { f -> "${label}:::${f}" }
                }
            )
            .mix(
                merge_nsforest_results_process.out.complete.map { meta, files ->
                    def label = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
                    def flist = files instanceof List ? files : [files]
                    flist.collect { f -> "${label}:::${f}" }
                }
            )
            .flatten()
            .collect()

        publish_results_process(
            params.organ,
            all_files_ch
        )
    } else {
        log.warn "WARNING: --github_token not set — skipping publish step"
    }
}
