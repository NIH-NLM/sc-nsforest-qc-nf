#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * main_plots_master.nf  -- STUB, not yet run
 *
 * Plots + master dataset summary only, for ONE dataset (Li, retina) whose
 * NSForest was run on the FULL data (Ajith's nsforest.pkl) and whose silhouette
 * was already computed on the full data. Reuses the modules of main.nf so
 * file names come out in the usual manner. main.nf is untouched.
 *
 * Skipped (their outputs come in as inputs): download, filter, dendrogram,
 * prep, run_nsforest, merge, compute_silhouette.
 *
 * Inputs (all paths may be S3 on Lifebit):
 *   --datasets_csv        the harvester CSV (use the Li row; first row is used)
 *   --organ               retina
 *   --custom_h5ad         the UNFILTERED h5ad Ajith's NSForest ran on (3,214,213 cells);
 *                         used for the NSForest plots and for the sha256 id
 *   --filtered_h5ad       the FILTERED h5ad the silhouette ran on (adata_filtered_..._cbd4aa);
 *                         used for the silhouette 2D projection
 *   --nsforest_results    results_ensg CSV  (TODO: made from nsforest.pkl, see below)
 *   --silhouette_scores, --cluster_summary, --annotation   existing silhouette outputs
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
include { viz_2D_projection_process }      from './modules/scsilhouette/viz_2D_projection.nf'
include { viz_distribution_process }       from './modules/scsilhouette/viz_distribution.nf'
include { viz_summary_process }            from './modules/scsilhouette/viz_summary.nf'
include { compute_summary_stats_process }  from './modules/scsilhouette/compute_summary_stats.nf'

params.datasets_csv       = null
params.organ              = null
params.custom_h5ad        = null
params.filtered_h5ad      = null
params.nsforest_results   = null
params.silhouette_scores  = null
params.cluster_summary    = null
params.annotation         = null
params.binary_scores      = null   // Ajith's supplementary-file.csv (binary scores, ENSG), passed through
params.outdir             = './'
params.publish_mode       = 'copy'
params.max_cells_per_cluster = 0
params.nsforest_seed         = 42

process checksum_h5ad_process {
    tag "sha256_${h5ad.name}"
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
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(binary)

    output:
    tuple val(meta), path("binary_scores_ensg_*.csv"), emit: binary_csv

    script:
    """
    cp ${binary} binary_scores_ensg_${meta.organ}_${meta.first_author}_${meta.journal.replace(' ', '_')}_${meta.year}_${meta.author_cell_type.replace(' ', '_')}_${meta.embedding.replace(' ', '_')}_${meta.dataset_version_id[-6..-1]}.csv
    """
}

workflow {
    ['datasets_csv', 'organ', 'custom_h5ad', 'filtered_h5ad', 'nsforest_results',
     'silhouette_scores', 'cluster_summary', 'annotation', 'binary_scores'].each { p ->
        if (!params[p]) { log.error "ERROR: --${p} is required"; exit 1 }
    }

    h5ad_ch = Channel.value(file(params.custom_h5ad))
    sha_out = checksum_h5ad_process(h5ad_ch)

    row_ch = Channel
        .fromPath(params.datasets_csv)
        .splitCsv(header: true, sep: ',')
        .first()

    // meta: same keys as main.nf; dataset_version_id replaced by the sha256 id
    meta_ch = row_ch.combine(sha_out.sha256).map { row, sha ->
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
            silhouette_input:                   'filtered',
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

    plots_process(h5ad_meta.join(results_ch))

    sil_ch = meta_ch.map { m ->
        tuple(m, file(params.silhouette_scores), file(params.cluster_summary), file(params.annotation))
    }

    viz_summary_process(
        sil_ch.join(results_ch).map { m, s, c, a, r -> tuple(m, s, c, a, r) }
    )
    viz_distribution_process(sil_ch)
    viz_2D_projection_process(filtered_meta)

    compute_summary_stats_process(
        sil_ch.join(results_ch).map { m, s, c, a, r -> tuple(m, s, c, a, r) }
    )
}
