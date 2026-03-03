#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { download_h5ad_process }                  from './modules/nsforest/download_h5ad.nf'
include { filter_adata_process }                   from './modules/nsforest/filter_adata.nf'
include { dendrogram_process }                     from './modules/nsforest/dendrogram.nf'
include { cluster_stats_process }                  from './modules/nsforest/cluster_stats.nf'
include { prep_medians_process }                   from './modules/nsforest/prep_medians.nf'
include { merge_medians_process }                  from './modules/nsforest/merge_medians.nf'
include { prep_binary_scores_process }             from './modules/nsforest/prep_binary_scores.nf'
include { merge_binary_scores_process }            from './modules/nsforest/merge_binary_scores.nf'
include { plot_histograms_process }                from './modules/nsforest/plot_histograms.nf'
include { run_nsforest_process }                   from './modules/nsforest/run_nsforest.nf'
include { merge_nsforest_results_process }         from './modules/nsforest/merge_nsforest_results.nf'
include { plots_process }                          from './modules/nsforest/plots.nf'
include { compute_silhouette_process }             from './modules/scsilhouette/compute_silhouette.nf'
include { viz_summary_process }                    from './modules/scsilhouette/viz_summary.nf'
include { viz_dotplot_process }                    from './modules/scsilhouette/viz_dotplot.nf'
include { viz_distribution_process }               from './modules/scsilhouette/viz_distribution.nf'
include { compute_summary_stats_process }          from './modules/scsilhouette/compute_summary_stats.nf'
// include { publish_results_process }             from './modules/publish/publish_results.nf'

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
                filter:           row.filter_normal
            ]
            tuple(meta, row.h5ad_url)
        }

    // Step 0a: Download h5ad from CellxGene URL
    downloaded_ch = download_h5ad_process(csv_rows_ch)

    // Step 0b: Filter
    filter_output_ch = filter_adata_process(
        downloaded_ch.h5ad
        .combine(uberon_ch)
        .combine(disease_ch)
        .combine(hsapdv_ch)

    )
    filter_output_ch = filter_adata_process(
        downloaded_ch.h5ad.combine(uberon_ch)
    )

    // Step 1: Dendrogram
    dendrogram_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )

    // Step 2: Cluster stats
    cluster_stats_output_ch = cluster_stats_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )

    // Step 3: Scatter by cluster
    scattered_clusters_ch = cluster_stats_output_ch.stats
        .flatMap { meta, h5ad, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, h5ad, row.cluster)
            }
        }

    // Step 4: Prep medians
    prep_medians_output_ch = prep_medians_process(scattered_clusters_ch)

    // Step 5: Gather medians
    merged_medians_ch = merge_medians_process(
        prep_medians_output_ch.partial.groupTuple()
    )

    // Step 6: Prep binary scores
    binary_scores_input_ch = merged_medians_ch.complete
        .map { meta, adata_prep, medians -> tuple(meta, adata_prep) }
        .combine(
            cluster_stats_output_ch.stats.map { meta, h5ad, stats -> tuple(meta, stats) },
            by: 0
        )
        .flatMap { meta, adata_prep, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }

    prep_binary_scores_output_ch = prep_binary_scores_process(binary_scores_input_ch)

    merged_binary_ch = merge_binary_scores_process(
        prep_binary_scores_output_ch.partial.groupTuple()
    )

    // Step 7: Plot histograms
    plot_histograms_process(
        merged_medians_ch.complete
            .map { meta, adata_prep, medians_files ->
                def medians_csv = medians_files instanceof List
                    ? medians_files.find { it.name.endsWith('.csv') } : medians_files
                tuple(meta, medians_csv)
            }
            .join(
                merged_binary_ch.complete.map { meta, adata_prep, binary_files ->
                    def binary_csv = binary_files instanceof List
                        ? binary_files.find { it.name.endsWith('.csv') } : binary_files
                    tuple(meta, binary_csv)
                }
            )
            .map { meta, medians_csv, binary_csv -> tuple(meta, medians_csv, binary_csv) }
    )

    // Step 8: Run NSForest
    nsforest_input_ch = merged_binary_ch.complete
        .map { meta, adata_prep, binary -> tuple(meta, adata_prep) }
        .combine(
            cluster_stats_output_ch.stats.map { meta, h5ad, stats -> tuple(meta, stats) },
            by: 0
        )
        .flatMap { meta, adata_prep, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }

    nsforest_output_ch = run_nsforest_process(nsforest_input_ch)

    merged_nsforest_ch = merge_nsforest_results_process(
        nsforest_output_ch.partial.groupTuple()
    )

    // Step 9: Plots
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

    // Step 10: Compute silhouette
    silhouette_output_ch = compute_silhouette_process(
        filter_output_ch.results
            .map { meta, h5ad, stats -> tuple(meta, h5ad) }
            .combine(uberon_ch)
    )

    // Step 11a: viz_summary
    viz_summary_process(
        silhouette_output_ch.results
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

    // Step 11b: viz_distribution
    viz_distribution_process(
        silhouette_output_ch.results.map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }
    )

    // Step 11c: viz_dotplot
    viz_dotplot_process(
        filter_output_ch.results.map { meta, h5ad, stats -> tuple(meta, h5ad) }
    )

    // Step 12: Summary statistics
    compute_summary_stats_process(
        silhouette_output_ch.results
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

//    publish_results_process(params.organ, nsforest_done_ch, silhouette_done_ch)
}
