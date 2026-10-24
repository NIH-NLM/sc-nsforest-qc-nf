#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * main_plots_master.nf  -- STUB, not yet run
 *
 * Plots + master dataset summary for ONE dataset (Li, retina) whose NSForest
 * was run on the FULL data (Ajith's nsforest.pkl). The silhouette scores are
 * computed here, over many tasks (modules/scsilhouette/silhouette_scatter.nf),
 * from the filtered h5ad. Reuses the modules of main.nf so file names come out
 * in the usual manner. main.nf is untouched.
 *
 * Skipped (their outputs come in as inputs): download, filter, dendrogram,
 * prep, run_nsforest, merge.
 *
 * Inputs (all paths may be S3 on Lifebit):
 *   --datasets_csv        the harvester CSV (use the Li row; first row is used)
 *   --organ               retina
 *   --custom_h5ad         the UNFILTERED h5ad Ajith's NSForest ran on (3,214,213 cells);
 *                         used for the NSForest plots and for the sha256 id
 *   --filtered_h5ad       the FILTERED h5ad the silhouette is computed on (adata_filtered_..._cbd4aa);
 *                         also used for the silhouette 2D projection
 *   --nsforest_results    results_ensg CSV  (TODO: made from nsforest.pkl, see below)
 *   --outdir              e.g. prod/data/retina/sc-nsforest-qc-nf/retina-li-2026-<id>
 *
 * Custom h5ad id: there is no cellxgene dataset_version_id for a custom file,
 * so the sha256 of the h5ad is computed and the FULL sha256 stands in for
 * dataset_version_id in meta. Like a cellxgene id, only its last 6 characters
 * appear in file names (get_output_prefix takes dataset_version_id[-6:]).
 * The full sha256 is written to <prefix>.h5ad.sha256 alongside the outputs.
 *
 * Open items (Anne decides):
 *   - binary scores: now exported from the pkl's AnnData varm (6946 genes x 140
 *     clusters, ENSG) as binary_scores_ensg.csv; pass that as --binary_scores
 *   - Ray's ETL needs: NSForest results (ENSG), master dataset summary, silhouette
 *     F-beta summary, binary scores (ENSG); exact names to confirm
 */

include { plots_process }                  from './modules/nsforest/plots.nf'
include { silhouette_scatter_gather }      from './modules/scsilhouette/silhouette_scatter.nf'
include { viz_2D_projection_process }      from './modules/scsilhouette/viz_2D_projection.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'
include { compute_summary_stats_process }  from './modules/scsilhouette/compute_summary_stats.nf'
include { plot_histograms_process }        from './modules/nsforest/plot_histograms.nf'
include { dendrogram_process }             from './modules/nsforest/dendrogram.nf'
include { cluster_stats_process }          from './modules/nsforest/cluster_stats.nf'
include { cluster_cid_mapping_process }    from './modules/nsforest/cluster_cid_mapping.nf'
include { publish_results_process }        from './modules/publish/publish_results.nf'

params.datasets_csv       = null
params.organ              = null
params.custom_h5ad        = null
params.filtered_h5ad      = null
params.nsforest_results   = null
params.binary_scores      = null   // binary scores (gene x cluster, ENSG) exported from Ajith's pkl, passed through
params.medians            = null   // medians (gene x cluster, ENSG) exported from Ajith's pkl, passed through
params.sha256             = null   // optional: the sha256 of --custom_h5ad, computed elsewhere; skips the checksum task
params.skip_plots         = false  // true: do not draw the NSForest plots (boxplots, scatter, dot, matrix, violin); they can come in a later run
params.outdir             = './'
params.github_token       = null   // without it nothing is published to GitHub
params.publish_repo       = 'NIH-NLM/nlm-ckn'
params.publish_dest_dir   = null   // default: data/prod/<organ>/sc-nsforest-qc-nf/results/<organ>-<author>-<journal>-<year>-<hash6>
params.publish_mode       = 'copy'
params.max_cells_per_cluster = 0
params.nsforest_seed         = 42

process checksum_h5ad_process {
    tag "sha256_${h5ad.name}"
    container 'ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:4c380184b7a60bbc7f26938b737130953a5b5982'
    // a large request, as for the dendrogram step, so that Batch picks a big machine with a fast network:
    // the 195 GB file is read in place, and a small machine reads it slowly
    cpus 8
    memory '384 GB'
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

// Pass-through: nothing is recomputed (no prep / NSForest / silhouette). The
// binary scores are renamed to the usual binary_scores_ensg_<prefix>.csv.
// Prefix follows get_output_prefix: organ_author_journal_year_clusterheader_embedding_vid6.
process passthrough_binary_scores_process {
    tag "binary_${meta.organ}_${meta.first_author}_${meta.year}_${meta.dataset_version_id}"
    container 'ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:4c380184b7a60bbc7f26938b737130953a5b5982'
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(binary)

    output:
    tuple val(meta), path("binary_scores_ensg_*.csv"), emit: binary_csv

    script:
    """
    cp ${binary} binary_scores_ensg_${meta.organ}_${meta.first_author}_${meta.journal.replace(' ', '_')}_${meta.year}_${meta.author_cell_type.replace(' ', '_')}_${meta.embedding.replace(' ', '_')}_${meta.dataset_version_id}.csv
    """
}

// Used when the sha256 of the unfiltered h5ad is given (--sha256): only writes the .sha256 file.
process write_sha256_process {
    tag "sha256_given"
    container 'ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:4c380184b7a60bbc7f26938b737130953a5b5982'
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

// Pass-through of Ajith's results (label 'results') and medians (label 'medians'):
// copied under the usual names results_ensg_<prefix>.csv and medians_ensg_<prefix>.csv.
process passthrough_ensg_process {
    tag "${label}_${meta.organ}_${meta.first_author}_${meta.year}_${meta.dataset_version_id}"
    container 'ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:4c380184b7a60bbc7f26938b737130953a5b5982'
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(source), val(label)

    output:
    tuple val(meta), path("*_ensg_*.csv"), emit: csv

    script:
    """
    cp ${source} ${label}_ensg_${meta.organ}_${meta.first_author}_${meta.journal.replace(' ', '_')}_${meta.year}_${meta.author_cell_type.replace(' ', '_')}_${meta.embedding.replace(' ', '_')}_${meta.dataset_version_id}.csv
    """
}

// master_s3_manifest.csv: every file this run publishes with its S3 address, plus one row
// for the filtered h5ad and one for the unfiltered h5ad, which this run only reads (the addresses given as
// --filtered_h5ad and --custom_h5ad). The h5ad files are not written to GitHub.
process manifest_process {
    tag "master_s3_manifest"
    container 'ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:4c380184b7a60bbc7f26938b737130953a5b5982'
    cpus 2
    memory '8 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    val file_names
    val s3_base
    val filtered_h5ad
    val unfiltered_h5ad

    output:
    path "master_s3_manifest.csv", emit: manifest

    script:
    def names_str = file_names.unique().sort().join('\n')
    def h5ad_name = filtered_h5ad.toString().tokenize('/').last()
    def unfiltered_name = unfiltered_h5ad.toString().tokenize('/').last()
    """
    cat > file_names.txt << 'NAMES_EOF'
${names_str}
NAMES_EOF
    nsforest-cli generate-s3-manifest --s3-base "${s3_base}"
    echo "${h5ad_name},${filtered_h5ad}" >> master_s3_manifest.csv
    echo "${unfiltered_name},${unfiltered_h5ad}" >> master_s3_manifest.csv
    """
}

workflow {
    ['datasets_csv', 'organ', 'custom_h5ad', 'filtered_h5ad', 'nsforest_results',
     'binary_scores', 'medians'].each { p ->
        if (!params[p]) { log.error "ERROR: --${p} is required"; exit 1 }
    }

    // the sha256 of the unfiltered h5ad is the dataset id: computed here, or given with --sha256
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

    row_ch = Channel
        .fromPath(params.datasets_csv)
        .splitCsv(header: true, sep: ',')
        .first()

    // meta: same keys as main.nf; dataset_version_id replaced by the sha256 id
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
            filtered_h5ad_path:                 params.filtered_h5ad,
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

    // TODO: process converting nsforest.pkl -> results_ensg CSV (+ markers.csv,
    // supplementary-file.csv) goes here once the pkl has been inspected.
    results_ch = meta_ch.map { m -> tuple(m, file(params.nsforest_results)) }
    h5ad_meta  = meta_ch.map { m -> tuple(m, file(params.custom_h5ad)) }
    filtered_meta = meta_ch.map { m -> tuple(m, file(params.filtered_h5ad)) }

    passthrough_binary_scores_process(meta_ch.map { m -> tuple(m, file(params.binary_scores)) })

    if (!params.skip_plots) {
        plots_process(h5ad_meta.join(results_ch))
    }

    // the silhouette scores, computed over many tasks from the filtered h5ad
    silhouette_output_ch = silhouette_scatter_gather(filtered_meta)

    sil_ch = silhouette_output_ch.scores
        .join(silhouette_output_ch.cluster_summary)
        .join(silhouette_output_ch.annotation)

    viz_summary_process(sil_ch.join(results_ch))
    viz_distribution_process(sil_ch)
    viz_2D_projection_process(filtered_meta)

    compute_summary_stats_process(sil_ch.join(results_ch))

    // Ajith's results and medians under the usual names
    passthrough_ensg_process(
        meta_ch.map { m -> tuple(m, file(params.nsforest_results), 'results') }
            .mix(meta_ch.map { m -> tuple(m, file(params.medians), 'medians') })
    )

    // the plots and tables that come from the medians, binary scores and the filtered file
    plot_histograms_process(meta_ch.map { m -> tuple(m, file(params.medians), file(params.binary_scores)) })
    dendrogram_process(filtered_meta)
    cluster_stats_process(filtered_meta)
    cluster_cid_mapping_process(filtered_meta)

    // the S3 manifest: the names of everything published above
    def s3_results_base = workflow.workDir.parent.toUriString() + '/results/' + params.outdir
    names_ch = Channel.empty()
        .mix(
            (params.skip_plots ? Channel.empty() : plots_process.out.plots),
            viz_summary_process.out.plots,
            viz_distribution_process.out.plots,
            viz_2D_projection_process.out.plots,
            compute_summary_stats_process.out.summary,
            passthrough_binary_scores_process.out.binary_csv,
            passthrough_ensg_process.out.csv,
            silhouette_output_ch.scores,
            silhouette_output_ch.cluster_summary,
            silhouette_output_ch.annotation,
            plot_histograms_process.out.histograms,
            dendrogram_process.out.svg,
            dendrogram_process.out.cluster_order,
            dendrogram_process.out.cluster_sizes,
            dendrogram_process.out.summary,
            cluster_stats_process.out.results,
            cluster_cid_mapping_process.out.results
        )
        .flatMap { meta, files -> (files instanceof List) ? files.flatten() : [files] }
        .map { f -> f.name }
        .mix(sha_file.map { f -> f.name })
        .collect()

    manifest_process(names_ch, s3_results_base, params.filtered_h5ad, params.custom_h5ad)

    // publish to a branch of the GitHub repository (only with --github_token)
    if (params.github_token) {
        publish_ch = Channel.empty()
            .mix(
                (params.skip_plots ? Channel.empty() : plots_process.out.plots),
                viz_summary_process.out.plots,
                viz_distribution_process.out.plots,
                viz_2D_projection_process.out.plots,
                compute_summary_stats_process.out.summary,
                passthrough_binary_scores_process.out.binary_csv,
                passthrough_ensg_process.out.csv,
                silhouette_output_ch.scores,
                silhouette_output_ch.cluster_summary,
                silhouette_output_ch.annotation,
                plot_histograms_process.out.histograms,
                dendrogram_process.out.svg,
                dendrogram_process.out.cluster_order,
                dendrogram_process.out.cluster_sizes,
                dendrogram_process.out.summary,
                cluster_stats_process.out.results,
                cluster_cid_mapping_process.out.results,
                meta_ch.combine(manifest_process.out.manifest)
            )
            .flatMap { meta, files -> ((files instanceof List) ? files.flatten() : [files]).collect { f -> tuple(meta, f) } }
            .groupTuple()
        publish_results_process(publish_ch)
    } else {
        log.warn "--github_token not set: the results are not published to GitHub"
    }
}
