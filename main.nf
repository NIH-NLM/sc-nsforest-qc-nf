#!/usr/bin/env nextflow

/*
 * sc-nsforest-qc-nf
 * =================
 * NSForest marker genes and silhouette scores for the datasets that
 * cellxgene-harvester-nf found, one run of the steps below for each dataset:
 *
 *   0a     download the h5ad (only when the CSV gives an https:// or s3:// address;
 *          a filtered h5ad written by cellxgene-harvester-nf is read as it is)
 *   0b     filter again by tissue, disease and age; drop clusters under min_cluster_size cells
 *   1      dendrogram, cluster statistics, cluster to cell ontology mapping
 *   2-3    medians and binary scores, histograms
 *   4-5    NSForest on batches of clusters, merge the results
 *   6      plots
 *   7-8    silhouette scores, summary and distribution plots, dataset summary
 *   9      the ETL JSON, the S3 manifest; publish to GitHub (needs --github_token)
 *
 * The input is the folder of <dataset_id>.filtered.json records from cellxgene-harvester-nf
 * (--harvester_json_dir); the curator adds author_cell_type and embedding to each record's
 * "curation" block. The output for the ETL is one JSON for each dataset: the whole record
 * plus the sc-nsforest-qc-nf results and the s3_ location of the final filtered h5ad.
 * All the choices are parameters; see nextflow.config and the README.
 */

include { cluster_cid_mapping_process }    from './modules/nsforest/cluster_cid_mapping.nf'
include { cluster_stats_process }          from './modules/nsforest/cluster_stats.nf'
include { compute_silhouette_process }     from './modules/scsilhouette/compute_silhouette.nf'
include { dendrogram_process }             from './modules/nsforest/dendrogram.nf'
include { download_h5ad_process }          from './modules/nsforest/download_h5ad.nf'
include { filter_adata_process }           from './modules/nsforest/filter_adata.nf'
include { generate_s3_manifest_process }   from './modules/publish/generate_s3_manifest.nf'
include { merge_nsforest_results_process } from './modules/nsforest/merge_nsforest_results.nf'
include { prep_process }                   from './modules/nsforest/prep.nf'
include { plot_histograms_process }        from './modules/nsforest/plot_histograms.nf'
include { plots_process }                  from './modules/nsforest/plots.nf'
include { publish_results_process }        from './modules/publish/publish_results.nf'
include { run_nsforest_process }           from './modules/nsforest/run_nsforest.nf'
include { compute_summary_stats_process }  from './modules/scsilhouette/compute_summary_stats.nf'
include { viz_2D_projection_process }      from './modules/scsilhouette/viz_2D_projection.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'
include { build_etl_json_process }         from './modules/publish/build_etl_json.nf'

// ---- the harvester records ------------------------------------------------------------

/* {'UBERON:1': 5, 'UBERON:2': 3} -> 'UBERON:1: 5; UBERON:2: 3' */
def id_summary(Map counts) {
    return counts.collect { id, n -> "${id}: ${n}" }.join('; ')
}

/*
 * The meta map of one dataset, read from its cellxgene-harvester-nf record
 * (<dataset_id>.filtered.json). The curator adds author_cell_type and embedding
 * to the record's "curation" block.
 */
def dataset_meta(rec) {
    def ds = rec.dataset
    return [
        organ:                      params.organ,
        dataset_id:                 ds.dataset_id,
        first_author:               ds.first_author.toString(),
        year:                       ds.year.toString(),
        author_cell_type:           rec.curation.author_cell_type,
        embedding:                  rec.curation.embedding,
        disease:                    rec.filtered_disease.join(' | '),
        filter_obs_column:          params.filter_obs_column,
        filter_obs_value:           params.filter_obs_value,
        doi:                        ds.doi,
        collection_name:            ds.collection_name,
        dataset_title:              ds.dataset_title,
        dataset_version_id:         ds.dataset_version_id,
        journal:                    ds.journal,
        collection_url:             ds.collection_url,
        explorer_url:               ds.explorer_url,
        h5ad_url:                   rec.filtered_h5ad_url,
        tissue_ontology_summary:    id_summary(rec.filtered_tissue_ontology_id_summary),
        assay_ontology_summary:     id_summary(rec.filtered_assay_ontology_id_summary),
        cell_type_ontology_summary: id_summary(rec.filtered_cell_type_ontology_id_summary),
        disease_ontology_summary:   id_summary(rec.filtered_disease_ontology_id_summary),
        sex_ontology_summary:       id_summary(rec.filtered_sex_ontology_id_summary),
        development_stage_summary:  id_summary(rec.filtered_development_stage_ontology_id_summary),
        session_id:                 workflow.sessionId.toString()[-6..-1],
    ]
}

// ---- the workflow ----------------------------------------------------------------------

workflow {

    // ---- check the parameters ------------------------------------------------------
    if (!params.harvester_json_dir) { error "Give --harvester_json_dir (the folder of <dataset_id>.filtered.json from cellxgene-harvester-nf)." }
    if (!params.organ)        { error "Give --organ." }
    if (!params.uberon_json)  { error "Give --uberon_json." }
    if (!params.disease_json) { error "Give --disease_json." }
    if (!params.hsapdv_json)  { error "Give --hsapdv_json." }

    log.info "sc-nsforest-qc-nf | organ ${params.organ} | records ${params.harvester_json_dir}"
    log.info "results in ${params.outdir} | work in ${workflow.workDir}"

    uberon_ch  = channel.value(file(params.uberon_json,  checkIfExists: true))
    disease_ch = channel.value(file(params.disease_json, checkIfExists: true))
    hsapdv_ch  = channel.value(file(params.hsapdv_json,  checkIfExists: true))

    // ---- the datasets -------------------------------------------------------------------
    def json_dir = file(params.harvester_json_dir, checkIfExists: true)

    csv_rows_ch = channel
        .fromPath("${json_dir}/*.filtered.json", checkIfExists: true)
        .map { json -> tuple(new groovy.json.JsonSlurper().parseText(json.text), json) }
        .filter { rec, _json ->
            def ref = rec.curation.reference?.toString()?.trim()?.toLowerCase()
            def who = "${rec.dataset.first_author} ${rec.dataset.year}"
            if (ref in ['exclude', 'delete', 'merge', 'question']) {
                log.info "Skipping ${who} — reference='${ref}'"
                return false
            }
            if (!(ref in ['yes', 'no', 'unk'])) {
                log.warn "Skipping ${who} — unrecognised reference value '${ref}'"
                return false
            }
            if (!rec.curation.author_cell_type?.toString()?.trim() || !rec.curation.embedding?.toString()?.trim()) {
                error "Dataset ${rec.dataset.dataset_id} (${who}) has an empty author_cell_type or embedding in its record; fill both in the \"curation\" block before running."
            }
            return true
        }
        .map { rec, json -> tuple(dataset_meta(rec), rec.filtered_h5ad_url, json) }

    // the JSON records, keyed by the same meta as every other output
    harvester_json_ch = csv_rows_ch
        .map { meta, _url, json -> tuple(meta, json) }

    // an h5ad address (https:// or s3://) is downloaded; anything else is a filtered h5ad
    // written by cellxgene-harvester-nf, found by file name in --h5ad_dir or next to the records
    h5ad_src_ch = csv_rows_ch
        .map { meta, url, _json -> tuple(meta, url) }
        .branch { _meta, url ->
            remote: url ==~ /^(https?|s3):\/\/.*/
            local:  true
        }

    local_h5ad_ch = h5ad_src_ch.local.map { meta, url ->
        def f = params.h5ad_dir ? file("${params.h5ad_dir}/${file(url).name}") : file("${json_dir}/${url}")
        if (!f.exists()) { error "Cannot find the filtered h5ad for ${meta.dataset_id}: ${f} (see --h5ad_dir)" }
        tuple(meta, f)
    }

    // Step 0a: Download h5ad from CellxGene URL
    downloaded_ch = download_h5ad_process(h5ad_src_ch.remote)
    source_h5ad_ch = downloaded_ch.h5ad.mix(local_h5ad_ch)

    // Step 0b: Filter — tissue + disease + age using per-row ontology term IDs
    filter_output_ch = filter_adata_process(
        source_h5ad_ch,
        uberon_ch,
        disease_ch,
        hsapdv_ch
    )

    // Convenience: filtered h5ad only channel
    filtered_h5ad_ch = filter_output_ch.h5ad

    // Step 1: Dendrogram
    dendrogram_output_ch = dendrogram_process(filtered_h5ad_ch)

    // Step 1b: Cluster statistics
    cluster_stats_process(filtered_h5ad_ch)

    // Step 1c: Cluster -> cell ontology ID mapping (4-column manual curation sheet)
    cluster_cid_mapping_process(filtered_h5ad_ch)

    // Step 2: Prep — medians + binary scores in ONE pass (single load / densification)
    prep_output_ch = prep_process(filtered_h5ad_ch)

    // Step 3: Plot histograms
    plot_histograms_process(
        prep_output_ch.medians_csv
            .join(prep_output_ch.binary_csv)
    )

    // Step 4: Scatter run_nsforest by cluster batch
    def batchSize = params.batch_size ?: 5

    nsforest_input_ch = filtered_h5ad_ch
        .join(prep_output_ch.medians_csv)
        .join(prep_output_ch.binary_csv)
        .join(dendrogram_output_ch.cluster_order)
        .flatMap { meta, h5ad, medians_csv, binary_csv, cluster_order_csv ->
            def clusters = cluster_order_csv
                .splitCsv(header: true)
                .collect { row -> row.cluster_order }
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

    // Step 9: the ETL JSON — the whole harvester record, the results and the s3_ location of the final h5ad
    etl_input_ch = compute_summary_stats_process.out.summary
        .map { meta, csvs -> tuple(meta.findAll { k, _v -> k != 'filtered_h5ad_path' }, csvs) }
        .join(harvester_json_ch)
        .join(filtered_h5ad_ch.map { meta, h5ad -> tuple(meta, h5ad.name) })
        .map { meta, csvs, json, h5ad_name -> tuple(meta, csvs, json, h5ad_name) }

    build_etl_json_process(etl_input_ch, uberon_ch, disease_ch, hsapdv_ch)

    // Step 9a: Publish + S3 Manifest
    def s3_results_base = workflow.workDir.parent.toUriString() + '/results'

    publish_base_ch = channel
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
            prep_output_ch.medians_csv,
            prep_output_ch.medians_csv_symbols,
            prep_output_ch.medians_pkl,
            prep_output_ch.medians_pkl_symbols,
            prep_output_ch.binary_csv,
            prep_output_ch.binary_csv_symbols,
            prep_output_ch.binary_pkl,
            prep_output_ch.binary_pkl_symbols,
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
            harvester_json_ch,
            build_etl_json_process.out.json,
        )
        .flatMap { meta, files ->
            def fileList = (files instanceof List) ? files.flatten() : [files]
            fileList.collect { f -> tuple(meta, f) }
        }

    // Step 9a: S3 manifest — collect names as strings, no file staging
    generate_s3_manifest_process(
        publish_base_ch
            .mix(filtered_h5ad_ch)
            .map { _meta, f -> f.name }
            .collect(),
        s3_results_base
    )

    // Step 9b: GitHub publish — conditional
    if (params.github_token) {
        publish_results_process(
            publish_base_ch
                .map { meta, f ->
                    def clean = meta.findAll { k, _v -> k != 'filtered_h5ad_path' }
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
