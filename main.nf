#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cluster_cid_mapping_process }    from './modules/nsforest/cluster_cid_mapping.nf'
include { cluster_stats_process }          from './modules/nsforest/cluster_stats.nf'
include { compute_silhouette_process }     from './modules/scsilhouette/compute_silhouette.nf'
include { dendrogram_process }             from './modules/nsforest/dendrogram.nf'
include { download_h5ad_process }          from './modules/nsforest/download_h5ad.nf'
include { filter_adata_process }           from './modules/nsforest/filter_adata.nf'
include { generate_s3_manifest_process }   from './modules/publish/generate_s3_manifest.nf'
include { merge_nsforest_results_process } from './modules/nsforest/merge_nsforest_results.nf'
include { prep_medians_process }           from './modules/nsforest/prep_medians.nf'
include { prep_binary_scores_process }     from './modules/nsforest/prep_binary_scores.nf'
include { plot_histograms_process }        from './modules/nsforest/plot_histograms.nf'
include { plots_process }                  from './modules/nsforest/plots.nf'
include { publish_results_process }        from './modules/publish/publish_results.nf'
include { run_nsforest_process }           from './modules/nsforest/run_nsforest.nf'
include { compute_summary_stats_process }  from './modules/scsilhouette/compute_summary_stats.nf'
include { viz_2D_projection_process }      from './modules/scsilhouette/viz_2D_projection.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'

params.batch_size        = 10
params.datasets_csv      = null
params.filter_obs_column = ''
params.filter_obs_value  = ''
params.organ             = null
params.uberon_json       = null
params.disease_json      = null
params.hsapdv_json       = null
params.min_cluster_size  = 5
params.outdir            = './'
params.publish_mode      = 'copy'

workflow {
    log.info "workflow.workDir    : ${workflow.workDir}"
    log.info "workflow.launchDir  : ${workflow.launchDir}"
    log.info "workflow.projectDir : ${workflow.projectDir}"
    log.info "workflow.runName    : ${workflow.runName}"
    log.info "params.outdir       : ${params.outdir}"

    if (!params.datasets_csv) { log.error "ERROR: --datasets_csv is required";  exit 1 }
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
                filter_obs_column:                  params.filter_obs_column,
                filter_obs_value:                   params.filter_obs_value,
                doi:                                row.doi,
                collection_name:                    row.collection_name,
                dataset_title:                      row.dataset_title,
                dataset_version_id:                 row.dataset_version_id,
                journal:                            row.journal,
                collection_url:                     row.collection_url,
                explorer_url:                       row.explorer_url,
                h5ad_url:                           row.h5ad_url,
                tissue_ontology_term_id:            row.tissue_ontology_term_id,
                disease_ontology_term_id:           row.disease_ontology_term_id,
                development_stage_ontology_term_id: row.development_stage_ontology_term_id,
                tissue_ontology_summary:            row.tissue_ontology_summary,
                assay_ontology_summary:             row.assay_ontology_summary,
                cell_type_ontology_summary:         row.cell_type_ontology_summary,
                disease_ontology_summary:           row.disease_ontology_summary,
                sex_ontology_summary:               row.sex_ontology_summary,
                development_stage_summary:          row.development_stage_summary,
                session_id:                         workflow.sessionId.toString()[-6..-1],
            ]
            tuple(meta, row.h5ad_url)
        }

    // Step 0a: Download h5ad from CellxGene URL
    downloaded_ch = download_h5ad_process(csv_rows_ch)

    // Step 0b: Filter — tissue + disease + age using per-row ontology term IDs
    filter_output_ch = filter_adata_process(
        downloaded_ch.h5ad,
        uberon_ch,
        disease_ch,
        hsapdv_ch
    )

    // Convenience: filtered h5ad only channel
    filtered_h5ad_ch = filter_output_ch.h5ad

    // Step 1: Dendrogram
    dendrogram_output_ch = dendrogram_process(filtered_h5ad_ch)

    // Step 1b: Cluster statistics
    cluster_stats_output_ch = cluster_stats_process(filtered_h5ad_ch)

    // Step 1c: Cluster -> cell ontology ID mapping (4-column manual curation sheet)
    cluster_cid_mapping_output_ch = cluster_cid_mapping_process(filtered_h5ad_ch)
    
    // Step 2a: Prep medians — runs once per dataset on full filtered h5ad
    prep_medians_output_ch = prep_medians_process(filtered_h5ad_ch)

    // Step 2b: Prep binary scores — runs once per dataset on full filtered h5ad
    prep_binary_scores_output_ch = prep_binary_scores_process(filtered_h5ad_ch)

    // Step 3: Plot histograms
    plot_histograms_process(
        prep_medians_output_ch.csv
            .join(prep_binary_scores_output_ch.csv)
    )

    // Step 4: Scatter run_nsforest by cluster batch
    def batchSize = params.batch_size ?: 10

    nsforest_input_ch = filtered_h5ad_ch
        .join(prep_medians_output_ch.csv)
        .join(prep_binary_scores_output_ch.csv)
        .join(dendrogram_output_ch.cluster_order)
        .flatMap { meta, h5ad, medians_csv, binary_csv, cluster_order_csv ->
            def clusters = cluster_order_csv
                .splitCsv(header: true)
                .collect { it.cluster_order }
            clusters.collate(batchSize).collect { batch ->
                tuple(meta, h5ad, medians_csv, binary_csv, batch.join(','))
            }
        }

    nsforest_output_ch = run_nsforest_process(nsforest_input_ch)

    // Step 5: Merge NSForest results (ENSG merge + symbol derivation from filtered h5ad)
    merge_input_ch = nsforest_output_ch.partial.groupTuple()
        .map { meta, file_lists -> tuple(meta, file_lists.flatten()) }
        .join(filtered_h5ad_ch)

    merged_nsforest_ch = merge_nsforest_results_process(merge_input_ch)
    
    // Step 6: Plots
    plots_process(
        filtered_h5ad_ch
            .join(merged_nsforest_ch.results_csv)
    )

    // Step 7: Compute silhouette
    silhouette_output_ch = compute_silhouette_process(filtered_h5ad_ch)

    // Step 8a: viz_summary
    viz_summary_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(merged_nsforest_ch.results_csv)
            .map { meta, scores, summary, annotation, nsforest_csv ->
                tuple(meta, scores, summary, annotation, nsforest_csv ?: file('NO_FILE'))
            }
    )
    
    // Step 8b: viz_distribution
    viz_distribution_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
    )

    // Step 8c: viz_2D_projection
    viz_2D_projection_process(filtered_h5ad_ch)

    // Step 8d: compute_summary_stats
    compute_summary_stats_process(
        filtered_h5ad_ch
            .join(silhouette_output_ch.scores)
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(merged_nsforest_ch.results_csv)
            .map { meta, h5ad, scores, cluster_summary, annotation, nsforest_csv ->
                def new_meta = meta + [filtered_h5ad_path: h5ad.toUriString()]
                tuple(new_meta, scores, cluster_summary, annotation, nsforest_csv ?: file('NO_FILE'))
            }
    )

    // Step 9: Publish + S3 Manifest
    def s3_results_base = workflow.workDir.parent.toUriString() + '/results'

    publish_base_ch = Channel
        .empty()
        .mix(
            dendrogram_process.out.cluster_order,
            dendrogram_process.out.cluster_sizes,
            dendrogram_process.out.summary,
            dendrogram_process.out.svg,
            cluster_stats_process.out.results,
            cluster_cid_mapping_process.out.results,
            filter_adata_process.out.cluster_sizes,
            filter_adata_process.out.cluster_order,
            filter_adata_process.out.summary,
            filter_adata_process.out.svg,
            plots_process.out.plots,
            prep_binary_scores_process.out.csv,
            prep_binary_scores_process.out.csv_symbols,
            prep_binary_scores_process.out.pkl,
            prep_binary_scores_process.out.pkl_symbols,
            prep_medians_process.out.csv,
            prep_medians_process.out.csv_symbols,
            prep_medians_process.out.pkl,
            prep_medians_process.out.pkl_symbols,
            merge_nsforest_results_process.out.results_csv,
            merge_nsforest_results_process.out.results_csv_symbols,
            merge_nsforest_results_process.out.results_pkl,
            merge_nsforest_results_process.out.results_pkl_symbols,
            merge_nsforest_results_process.out.markers,
            merge_nsforest_results_process.out.markers_symbols,
            merge_nsforest_results_process.out.markers_ontarget,
            merge_nsforest_results_process.out.markers_ontarget_symbols,
            merge_nsforest_results_process.out.markers_ontarget_supp,
            merge_nsforest_results_process.out.markers_ontarget_supp_symbols,
            merge_nsforest_results_process.out.gene_selection,
            merge_nsforest_results_process.out.gene_selection_symbols,
            plot_histograms_process.out.histograms,
            compute_silhouette_process.out.scores,
            compute_silhouette_process.out.cluster_summary,
            compute_silhouette_process.out.annotation,
            viz_2D_projection_process.out.plots,
            viz_distribution_process.out.plots,
            viz_summary_process.out.plots,
            compute_summary_stats_process.out.summary,
        )
        .flatMap { meta, files ->
            def fileList = (files instanceof List) ? files.flatten() : [files]
            fileList.collect { f -> tuple(meta, f) }
        }

    // Step 9a: S3 manifest — always, includes h5ad, runs once across all datasets
    generate_s3_manifest_process(
        publish_base_ch
            .mix(filtered_h5ad_ch)
            .map { meta, f -> f }
            .collect(),
        s3_results_base
    )

    // Step 9b: GitHub publish — conditional
    if (params.github_token) {
        publish_results_process(
            publish_base_ch
                .map { meta, f ->
                    def clean = meta.findAll { k, v -> k != 'filtered_h5ad_path' }
                    tuple(clean, f)
                }
                .groupTuple()
                .map { meta, file_lists -> tuple(meta, file_lists.flatten()) }
                .combine(generate_s3_manifest_process.out.manifest)
                .map { meta, files, manifest -> tuple(meta, files + [manifest]) }
        )
    } else {
        log.warn "WARNING: --github_token not set -- skipping publish step"
    }
}
