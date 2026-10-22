#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * main.nf (branch neocortex-run) -- NS-Forest markers, marker combinations and silhouette scores
 * for ONE custom h5ad with several cluster headers: the Jorstad et al. 2023 human neocortex atlas,
 * downsampled to 73,144 cells (Zenodo 10.5281/zenodo.13172641).
 *
 * The h5ad is used as it is: no download, no filter. One run of the steps below for each row of the
 * datasets CSV (one row per cluster header: CrossArea_subclass and CrossArea_cluster):
 *
 *   1      dendrogram (cluster order), cluster statistics, cluster to cell ontology mapping
 *   2-3    medians and binary scores, histograms
 *   4-5    NS-Forest on batches of clusters, merge the results
 *   6      marker combinations (DecisionTree, combinations=True): the cluster header, with the
 *          NS-Forest markers of its parent subclasses (--subclass_header)
 *   7      plots
 *   8      silhouette scores as scatter-gather (see modules/scsilhouette/silhouette_scatter.nf),
 *          summary and distribution plots, dataset summary
 *   9      the S3 manifest; publish to GitHub (needs --github_token)
 *
 * No ETL JSON in this run: the CSV files and plots are what the ETL reads.
 *
 * Inputs (paths may be S3 on Lifebit):
 *   --datasets_csv   one row per cluster header (see params/neocortex_datasets.csv)
 *   --organ          neocortex
 *   --custom_h5ad    the downsampled h5ad
 *   --outdir         e.g. prod/data/neocortex/sc-nsforest-qc-nf/neocortex-jorstad-2023
 *
 * Dataset id: a custom file has no cellxgene dataset_version_id, so the sha256 of the h5ad stands in
 * for it. Like a cellxgene id, only its last 6 characters appear in file names. The full sha256 is
 * written to <h5ad name>.sha256 beside the outputs.
 */

include { cluster_cid_mapping_process }    from './modules/nsforest/cluster_cid_mapping.nf'
include { cluster_stats_process }          from './modules/nsforest/cluster_stats.nf'
include { combinations_process }           from './modules/nsforest/combinations.nf'
include { compare_published_process as compare_subclass_process }      from './modules/nsforest/compare_published.nf'
include { compare_published_process as compare_combinations_process }  from './modules/nsforest/compare_published.nf'
include { compute_summary_stats_process }  from './modules/scsilhouette/compute_summary_stats.nf'
include { dendrogram_process }             from './modules/nsforest/dendrogram.nf'
include { merge_nsforest_results_process } from './modules/nsforest/merge_nsforest_results.nf'
include { plot_histograms_process }        from './modules/nsforest/plot_histograms.nf'
include { plots_process }                  from './modules/nsforest/plots.nf'
include { prep_process }                   from './modules/nsforest/prep.nf'
include { publish_results_process }        from './modules/publish/publish_results.nf'
include { run_nsforest_process }           from './modules/nsforest/run_nsforest.nf'
include { silhouette_scatter_gather }      from './modules/scsilhouette/silhouette_scatter.nf'
include { viz_2D_projection_process }      from './modules/scsilhouette/viz_2D_projection.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'

params.datasets_csv       = null
params.organ              = null
params.custom_h5ad        = null
params.subclass_header    = 'CrossArea_subclass'   // obs column of the parent subclasses
params.combinations_header = 'CrossArea_cluster'   // the cluster header whose combinations are evaluated
params.published_subclass_csv     = null   // optional: the published subclass results (Zenodo); compared with our subclass results
params.published_combinations_csv = null   // optional: the published cluster combinatorial results; compared with ours
params.sha256             = null   // optional: the sha256 of --custom_h5ad, computed elsewhere; skips the checksum task
params.skip_plots         = false  // true: do not draw the NS-Forest plots (boxplots, scatter, dot, matrix, violin)
params.outdir             = './'
params.github_token       = null   // without it nothing is published to GitHub
params.publish_repo       = 'NIH-NLM/nlm-ckn'
params.publish_dest_dir   = null   // default: data/prod/<organ>/sc-nsforest-qc-nf/results/<organ>-<author>-<journal>-<year>-<hash6>
params.publish_mode       = 'copy'
params.max_cells_per_cluster = 0
params.nsforest_seed         = 42

process checksum_h5ad_process {
    tag "sha256_${h5ad.name}"
    container params.nsforest_image
    cpus 4
    memory '32 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path h5ad

    output:
    env  'SHA256',               emit: sha256
    path "${h5ad.name}.sha256",  emit: file

    script:
    """
    SHA256=\$(sha256sum ${h5ad} | cut -d' ' -f1)
    echo "\$SHA256  ${h5ad.name}" > ${h5ad.name}.sha256
    """
}

// Used when the sha256 of the h5ad is given (--sha256): only writes the .sha256 file.
process write_sha256_process {
    tag "sha256_given"
    container params.nsforest_image
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    val sha256
    val h5ad_name

    output:
    path "${h5ad_name}.sha256", emit: file

    script:
    """
    echo "${sha256}  ${h5ad_name}" > ${h5ad_name}.sha256
    """
}

// master_s3_manifest.csv: every file this run publishes with its S3 address, plus one row for the h5ad,
// which this run only reads (the address given as --custom_h5ad). The h5ad is not written to GitHub.
process manifest_process {
    tag "master_s3_manifest"
    container params.nsforest_image
    cpus 2
    memory '8 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    val file_names
    val s3_base
    val h5ad

    output:
    path "master_s3_manifest.csv", emit: manifest

    script:
    def names_str = file_names.unique().sort().join('\n')
    def h5ad_name = h5ad.toString().tokenize('/').last()
    """
    cat > file_names.txt << 'NAMES_EOF'
${names_str}
NAMES_EOF
    nsforest-cli generate-s3-manifest --s3-base "${s3_base}"
    echo "${h5ad_name},${h5ad}" >> master_s3_manifest.csv
    """
}

workflow {
    ['datasets_csv', 'organ', 'custom_h5ad'].each { p ->
        if (!params[p]) { log.error "ERROR: --${p} is required"; exit 1 }
    }

    // the sha256 of the h5ad is the dataset id: computed here, or given with --sha256
    if (params.sha256) {
        if (!(params.sha256.toString().trim().toLowerCase() ==~ /[0-9a-f]{64}/)) {
            log.error "ERROR: --sha256 must be 64 hexadecimal characters"; exit 1
        }
        sha_value = Channel.value(params.sha256.toString().trim().toLowerCase())
        sha_file  = write_sha256_process(sha_value, params.custom_h5ad.toString().tokenize('/').last()).file
    } else {
        sha_out   = checksum_h5ad_process(Channel.value(file(params.custom_h5ad)))
        sha_value = sha_out.sha256
        sha_file  = sha_out.file
    }

    // one row for each cluster header
    row_ch = Channel
        .fromPath(params.datasets_csv)
        .splitCsv(header: true, sep: ',')

    // meta: the keys the modules use; dataset_version_id is the sha256 id
    meta_ch = row_ch.combine(sha_value).map { row, sha ->
        [
            organ:                              params.organ,
            first_author:                       row.first_author,
            year:                               row.year,
            author_cell_type:                   row.author_cell_type,
            embedding:                          row.embedding,
            disease:                            row.disease,
            filter:                             row.filter_normal,
            filter_obs_column:                  '',
            filter_obs_value:                   '',
            doi:                                row.doi,
            collection_name:                    row.collection_name,
            dataset_title:                      row.dataset_title,
            dataset_version_id:                 sha,
            journal:                            row.journal,
            collection_url:                     row.collection_url,
            explorer_url:                       row.explorer_url,
            h5ad_url:                           params.custom_h5ad,
            filtered_h5ad_path:                 params.custom_h5ad,
            nsforest_input:                     'unfiltered',
            silhouette_input:                   'unfiltered',
            tissue_ontology_term_id:            row.tissue_ontology_term_id,
            disease_ontology_term_id:           row.disease_ontology_term_id,
            development_stage_ontology_term_id: row.development_stage_ontology_term_id,
            tissue_ontology_summary:            row.tissue_ontology_summary,
            assay_ontology_summary:             row.assay_ontology_summary,
            cell_type_ontology_summary:         row.cell_type_ontology_summary,
            disease_ontology_summary:           row.disease_ontology_summary,
            sex_ontology_summary:               row.sex_ontology_summary,
            development_stage_summary:          row.development_stage_summary,
            sha256:                             sha,
            session_id:                         workflow.sessionId.toString()[-6..-1],
        ]
    }

    h5ad_ch = meta_ch.map { m -> tuple(m, file(params.custom_h5ad)) }

    // Step 1: dendrogram (cluster order), cluster statistics, cluster -> cell ontology mapping
    dendrogram_output_ch = dendrogram_process(h5ad_ch)
    cluster_stats_process(h5ad_ch)
    cluster_cid_mapping_process(h5ad_ch)

    // Step 2: medians and binary scores in one pass
    prep_output_ch = prep_process(h5ad_ch)

    // Step 3: histograms
    plot_histograms_process(
        prep_output_ch.medians_csv
            .join(prep_output_ch.binary_csv)
    )

    // Step 4: NS-Forest on batches of clusters
    def batchSize = params.batch_size ?: 5

    nsforest_input_ch = h5ad_ch
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

    // Step 5: merge the batches
    merge_input_ch = nsforest_output_ch.partial.groupTuple()
        .map { meta, file_lists -> tuple(meta, file_lists.flatten()) }
        .join(h5ad_ch)

    merged_nsforest_ch = merge_nsforest_results_process(merge_input_ch)

    // Step 6: marker combinations for the cluster header, with the markers of its parent subclasses.
    // Both rows share the same dataset id, which pairs them.
    subclass_results_ch = merged_nsforest_ch.results_csv
        .filter { meta, f -> meta.author_cell_type == params.subclass_header }
        .map { meta, f -> tuple(meta.dataset_version_id, f) }

    combinations_input_ch = merged_nsforest_ch.results_csv
        .filter { meta, f -> meta.author_cell_type == params.combinations_header }
        .map { meta, f -> tuple(meta.dataset_version_id, meta, f) }
        .join(subclass_results_ch)
        .map { vid, meta, cluster_results, subclass_results ->
            tuple(meta, file(params.custom_h5ad), params.subclass_header, subclass_results, cluster_results)
        }

    combinations_process(combinations_input_ch)

    // Step 6b: compare with the published results, cluster by cluster (only when they are given)
    comparison_ch = Channel.empty()
    if (params.published_subclass_csv) {
        compare_subclass_process(
            merged_nsforest_ch.results_csv
                .filter { meta, f -> meta.author_cell_type == params.subclass_header }
                .map { meta, f -> tuple(meta, f, file(params.published_subclass_csv), 'published_subclass') }
        )
        comparison_ch = comparison_ch.mix(compare_subclass_process.out.comparison)
    }
    if (params.published_combinations_csv) {
        compare_combinations_process(
            combinations_process.out.results
                .map { meta, f -> tuple(meta, f, file(params.published_combinations_csv), 'published_combinations') }
        )
        comparison_ch = comparison_ch.mix(compare_combinations_process.out.comparison)
    }

    // Step 7: plots
    if (!params.skip_plots) {
        plots_process(h5ad_ch.join(merged_nsforest_ch.results_csv))
    }

    // Step 8: silhouette (scatter-gather) and its plots and summary
    silhouette_output_ch = silhouette_scatter_gather(h5ad_ch)

    viz_summary_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(merged_nsforest_ch.results_csv)
    )
    viz_distribution_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
    )
    viz_2D_projection_process(h5ad_ch)

    compute_summary_stats_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(merged_nsforest_ch.results_csv)
    )

    // Step 9: everything the run writes, for the S3 manifest and for GitHub
    def s3_results_base = workflow.workDir.parent.toUriString() + '/results/' + params.outdir

    outputs_ch = Channel.empty()
        .mix(
            dendrogram_output_ch.svg,
            dendrogram_output_ch.cluster_order,
            dendrogram_output_ch.cluster_sizes,
            dendrogram_output_ch.summary,
            cluster_stats_process.out.results,
            cluster_cid_mapping_process.out.results,
            prep_output_ch.medians_csv,
            prep_output_ch.medians_csv_symbols,
            prep_output_ch.medians_pkl,
            prep_output_ch.medians_pkl_symbols,
            prep_output_ch.binary_csv,
            prep_output_ch.binary_csv_symbols,
            prep_output_ch.binary_pkl,
            prep_output_ch.binary_pkl_symbols,
            plot_histograms_process.out.histograms,
            merged_nsforest_ch.results_csv,
            merged_nsforest_ch.results_csv_symbols,
            merged_nsforest_ch.results_pkl,
            merged_nsforest_ch.results_pkl_symbols,
            merged_nsforest_ch.markers,
            merged_nsforest_ch.markers_symbols,
            merged_nsforest_ch.markers_ontarget,
            merged_nsforest_ch.markers_ontarget_symbols,
            merged_nsforest_ch.markers_ontarget_supp,
            merged_nsforest_ch.markers_ontarget_supp_symbols,
            merged_nsforest_ch.gene_selection,
            merged_nsforest_ch.gene_selection_symbols,
            combinations_process.out.results,
            comparison_ch,
            (params.skip_plots ? Channel.empty() : plots_process.out.plots),
            silhouette_output_ch.scores,
            silhouette_output_ch.cluster_summary,
            silhouette_output_ch.annotation,
            viz_summary_process.out.plots,
            viz_distribution_process.out.plots,
            viz_2D_projection_process.out.plots,
            compute_summary_stats_process.out.summary
        )

    names_ch = outputs_ch
        .flatMap { meta, files -> (files instanceof List) ? files.flatten() : [files] }
        .map { f -> f.name }
        .mix(sha_file.map { f -> f.name })
        .collect()

    manifest_process(names_ch, s3_results_base, params.custom_h5ad)

    // publish to a branch of the GitHub repository (only with --github_token)
    if (params.github_token) {
        publish_ch = outputs_ch
            .mix(meta_ch.combine(manifest_process.out.manifest))
            .flatMap { meta, files -> ((files instanceof List) ? files.flatten() : [files]).collect { f -> tuple(meta, f) } }
            .groupTuple()
        publish_results_process(publish_ch)
    } else {
        log.warn "--github_token not set: the results are not published to GitHub"
    }
}
