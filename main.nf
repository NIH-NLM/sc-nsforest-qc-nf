#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// NSForest processes
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

// scsilhouette processes
include { compute_silhouette_process }             from './modules/scsilhouette/compute_scsilhouette.nf'
include { viz_summary_process }                    from './modules/scsilhouette/viz_summary.nf'
include { viz_dotplot_process }                    from './modules/scsilhouette/viz_dotplot.nf'
include { viz_distribution_process }               from './modules/scsilhouette/viz_distribution.nf'
include { compute_summary_stats_process }          from './modules/scsilhouette/compute_summary_stats.nf'

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

params.datasets_csv  = null    // Path to homo_sapiens_{organ}_harvester_final.csv
params.organ         = null    // e.g. "heart" — used for output directory naming
params.uberon_json   = null    // Path to uberon_{organ}.json from cellxgene-harvester resolve-uberon
params.min_age       = 15      // Minimum donor age in years for adult cell filtering
params.min_cluster_size = 5    // Minimum cells per cluster — smaller clusters dropped and logged
params.outdir        = './results'
params.publish_mode  = 'copy'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // -----------------------------------------------------------------------
    // Validate required parameters
    // -----------------------------------------------------------------------

    if (!params.datasets_csv) {
        log.error "ERROR: --datasets_csv is required (homo_sapiens_{organ}_harvester_final.csv)"
        exit 1
    }

    if (!params.organ) {
        log.error "ERROR: --organ is required (e.g. heart, kidney)"
        exit 1
    }

    if (!params.uberon_json) {
        log.error "ERROR: --uberon_json is required (uberon_{organ}.json from cellxgene-harvester resolve-uberon)"
        exit 1
    }

    // -----------------------------------------------------------------------
    // Stage the uberon JSON as a Nextflow file channel.
    // One file, shared across all parallel h5ad jobs.
    // -----------------------------------------------------------------------

    uberon_ch = Channel.fromPath(params.uberon_json, checkIfExists: true)

    // -----------------------------------------------------------------------
    // Read the harvester CSV and build the meta map for each dataset row.
    //
    // CSV columns expected (from homo_sapiens_{organ}_harvester_final.csv):
    //   h5ad_file, first_author, year, author_cell_type, embedding,
    //   disease, filter_normal
    //
    // Note: tissue filtering is now driven entirely by uberon_json (which
    // encodes the full anatomical unit — e.g. heart + pericardium + all
    // descendant UBERON terms). params.tissue is no longer needed.
    // -----------------------------------------------------------------------

    csv_rows_ch = Channel
        .fromPath(params.datasets_csv)
        .ifEmpty { exit 1, "Cannot find datasets CSV: ${params.datasets_csv}" }
        .splitCsv(header: true, sep: ',')
        .filter { row ->
            def ref = row.reference?.trim().toLowerCase()
            if (ref in ['question', 'merge']) {
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
                organ:           params.organ,
                first_author:    row.first_author,
                year:            row.year,
                author_cell_type: row.author_cell_type,
                embedding:       row.embedding,
                disease:         row.disease,
                filter:          row.filter_normal
            ]
            tuple(meta, file(row.h5ad_file))
        }

    // -----------------------------------------------------------------------
    // Step 0: Filter — tissue (UBERON IDs) + disease (PATO:0000461) + age
    //
    // Each h5ad is combined with the single uberon_json file so Nextflow
    // stages it into every parallel filter job's work directory.
    // -----------------------------------------------------------------------

    filter_input_ch = csv_rows_ch.combine(uberon_ch)
    // Channel shape: (meta, h5ad, uberon_json)

    filter_output_ch = filter_adata_process(filter_input_ch)
    // Emits: (meta, adata_filtered.h5ad, [stats/svg files])

    // Extract just (meta, filtered_h5ad) for downstream processes
    filtered_h5ad_ch = filter_output_ch.results
        .map { meta, filtered_h5ad, stats -> tuple(meta, filtered_h5ad) }

    // -----------------------------------------------------------------------
    // Step 1: Dendrogram (on filtered cells)
    // -----------------------------------------------------------------------

    dendrogram_output_ch = dendrogram_process(filtered_h5ad_ch)

    // -----------------------------------------------------------------------
    // Step 2: Cluster stats — drives scatter/gather parallelisation
    // -----------------------------------------------------------------------

    cluster_stats_output_ch = cluster_stats_process(filtered_h5ad_ch)

    // -----------------------------------------------------------------------
    // Step 3: SCATTER — one channel item per cluster
    // -----------------------------------------------------------------------

    scattered_clusters_ch = cluster_stats_output_ch
        .flatMap { meta, h5ad, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, h5ad, row.cluster)
            }
        }

    // -----------------------------------------------------------------------
    // Step 4: Prep medians — parallel by cluster
    // -----------------------------------------------------------------------

    prep_medians_output_ch = prep_medians_process(scattered_clusters_ch)

    // -----------------------------------------------------------------------
    // Step 5: GATHER medians — group all partial files back by dataset
    // -----------------------------------------------------------------------

    gathered_medians_ch = prep_medians_output_ch.groupTuple()

    merged_medians_ch = merge_medians_process(gathered_medians_ch)
    // Emits: (meta, adata_prep.h5ad, medians.{csv,pkl})

    // -----------------------------------------------------------------------
    // Step 6: Prep binary scores — parallel by cluster (reuse cluster list)
    // -----------------------------------------------------------------------

    binary_scores_input_ch = merged_medians_ch.complete
        .map { meta, adata_prep, medians -> tuple(meta, adata_prep) }
        .combine(cluster_stats_output_ch.map { meta, h5ad, stats -> stats }, by: 0)
        .flatMap { meta, adata_prep, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }

    prep_binary_scores_output_ch = prep_binary_scores_process(binary_scores_input_ch)

    gathered_binary_ch = prep_binary_scores_output_ch.groupTuple()

    merged_binary_ch = merge_binary_scores_process(gathered_binary_ch)
    // Emits: (meta, adata_prep.h5ad, binary_scores.{csv,pkl})

    // -----------------------------------------------------------------------
    // Step 7: Plot histograms
    // -----------------------------------------------------------------------

    histograms_input_ch = merged_medians_ch.complete
        .map { meta, adata_prep, medians_files ->
            // medians_files is a list [medians.csv, medians.pkl] — take the csv
            def medians_csv = medians_files instanceof List
                ? medians_files.find { it.name.endsWith('.csv') }
                : medians_files
            tuple(meta, medians_csv)
        }
        .join(
            merged_binary_ch.complete.map { meta, adata_prep, binary_files ->
                def binary_csv = binary_files instanceof List
                    ? binary_files.find { it.name.endsWith('.csv') }
                    : binary_files
                tuple(meta, binary_csv)
            }
        )
        .map { meta, medians_csv, binary_csv -> tuple(meta, medians_csv, binary_csv) }

    plot_histograms_process(histograms_input_ch)

    // -----------------------------------------------------------------------
    // Step 8: Run NSForest — parallel by cluster
    // -----------------------------------------------------------------------

    nsforest_input_ch = merged_binary_ch.complete
        .map { meta, adata_prep, binary -> tuple(meta, adata_prep) }
        .combine(cluster_stats_output_ch.map { meta, h5ad, stats -> stats }, by: 0)
        .flatMap { meta, adata_prep, stats_csv ->
            stats_csv.splitCsv(header: true).collect { row ->
                tuple(meta, adata_prep, row.cluster)
            }
        }

    nsforest_output_ch = run_nsforest_process(nsforest_input_ch)

    gathered_nsforest_ch = nsforest_output_ch.groupTuple()

    merged_nsforest_ch = merge_nsforest_results_process(gathered_nsforest_ch)
    // Emits: (meta, results.{csv,pkl})

    // -----------------------------------------------------------------------
    // Step 9: Plots (boxplots, scatter, expression)
    // -----------------------------------------------------------------------

    plots_input_ch = filtered_h5ad_ch
        .join(merged_nsforest_ch.complete.map { meta, results_files ->
            def results_csv = results_files instanceof List
                ? results_files.find { it.name.endsWith('.csv') }
                : results_files
            tuple(meta, results_csv)
        })
        .map { meta, h5ad, results_csv -> tuple(meta, h5ad, results_csv) }

    plots_process(plots_input_ch)

    // -----------------------------------------------------------------------
    // Step 10: Compute silhouette scores
    //
    // Uses the filtered h5ad (already tissue/disease/age filtered by step 0).
    // The uberon_json and min_age are passed through so the python code can
    // record filter parameters in annotation.json for provenance.
    // filter-normal flag re-applies the same filter on top of the already
    // filtered h5ad — this is intentional for datasets where filter_adata
    // was skipped or to guarantee consistency.
    // -----------------------------------------------------------------------

    silhouette_input_ch = filtered_h5ad_ch.combine(uberon_ch)
    // Channel shape: (meta, adata_filtered.h5ad, uberon_json)

    silhouette_output_ch = compute_silhouette_process(silhouette_input_ch)
    // Emits: (meta, silhouette_scores.csv, cluster_summary.csv, annotation.json)

    // -----------------------------------------------------------------------
    // Step 11: Silhouette visualizations
    // -----------------------------------------------------------------------

    // viz_summary — with optional NSForest F-scores
    viz_summary_input_ch = silhouette_output_ch.results
        .join(
            merged_nsforest_ch.complete.map { meta, results_files ->
                def results_csv = results_files instanceof List
                    ? results_files.find { it.name.endsWith('.csv') }
                    : results_files
                tuple(meta, results_csv)
            },
            remainder: true
        )
        .map { meta, scores, summary, annotation, nsforest_csv ->
            tuple(meta, scores, summary, annotation,
                  nsforest_csv ?: file('NO_FILE'))
        }

    viz_summary_process(viz_summary_input_ch)

    // viz_distribution
    viz_distribution_input_ch = silhouette_output_ch.results
        .map { meta, scores, summary, annotation ->
            tuple(meta, scores, summary, annotation)
        }

    viz_distribution_process(viz_distribution_input_ch)

    // viz_dotplot (uses filtered h5ad)
    viz_dotplot_process(filtered_h5ad_ch)

    // -----------------------------------------------------------------------
    // Step 12: Dataset-level summary statistics
    // -----------------------------------------------------------------------

    summary_stats_input_ch = silhouette_output_ch.results
        .join(
            merged_nsforest_ch.complete.map { meta, results_files ->
                def results_csv = results_files instanceof List
                    ? results_files.find { it.name.endsWith('.csv') }
                    : results_files
                tuple(meta, results_csv)
            },
            remainder: true
        )
        .map { meta, scores, summary, annotation, nsforest_csv ->
            tuple(meta, scores, summary, annotation,
                  nsforest_csv ?: file('NO_FILE'))
        }

    compute_summary_stats_process(summary_stats_input_ch)
}
