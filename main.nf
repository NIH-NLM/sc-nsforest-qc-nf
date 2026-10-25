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
 *   --outdir         optional; the default is results, as in the other runs
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
include { silhouette_scatter_over_datasets } from './modules/scsilhouette/silhouette_scatter.nf'
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
params.outdir             = 'results'
params.github_token       = null   // without it nothing is published to GitHub
params.publish_repo       = 'NIH-NLM/nlm-ckn'
params.publish_dest_dir   = null   // default: data/prod/<organ>/sc-nsforest-qc-nf/results/<organ>-<author>-<journal>-<year>-<hash6>
params.publish_mode       = 'copy'
params.max_cells_per_cluster = 0       // for the NS-Forest runs (the calculation is never subsampled)
params.plots_max_cells_per_cluster = 0   // for the plots only; 0 = all the cells
params.nsforest_seed         = 42
params.silhouette_datasets_csv = null   // one row per supercluster h5ad: order,supercluster,dataset_version_id,h5ad_url,expected_cells
params.renee_results      = null   // the NS-Forest results table of the paper (clusterSize, sampled_clusterSize), passed through
params.dend_order         = null   // dendOrder_author_cell_type.csv (column CrossArea_cluster_order), passed through
params.publish_header     = 'CrossArea_cluster'   // the cluster header whose results are published (one copy of every file)

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

// The NS-Forest results of the paper (Renee's table: clusterSize, sampled_clusterSize, PPV, ...) are passed through, as
// Ajith's results are for retina; our own combinations result is only the CHECK that they are the same (the run stops
// if a marker, a score or a count differs). With them the ENSG and the gene-symbol versions of the results, the
// markers, the binary scores and the medians are made under the usual names (csv and pkl). The ENSG ids are the ones of
// the h5ad itself: raw.var_names has the Ensembl id of every gene, in the order of var_names (the symbols).
// The results get one column, software_version: the calculation was NS-Forest 4.0.0, as for the paper.
// plots_input_results_*.csv is the same table with PPV called precision, for the plots; it is not published.
process passthrough_renee_process {
    tag "results_${meta.organ}_${meta.first_author}_${meta.year}_${meta.dataset_version_id}"
    container params.nsforest_image
    cpus 2
    memory '32 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(renee, stageAs: 'in_renee.csv'), path(binary, stageAs: 'in_binary.csv'), path(medians, stageAs: 'in_medians.csv'), path(h5ad), path(combos, stageAs: 'in_combinations.csv')

    output:
    tuple val(meta), path("*_{ensg,symbols}_${meta.organ}_*.{csv,pkl}"), emit: files
    tuple val(meta), path("plots_input_results_*.csv"),                  emit: plots_results
    path "renee_vs_ours_*.log",                                          emit: check

    script:
    def vid    = meta.dataset_version_id.toString()[-6..-1]
    def prefix = "${meta.organ}_${meta.first_author}_${meta.journal.replace(' ', '_')}_${meta.year}_${meta.author_cell_type.replace(' ', '_')}_${meta.embedding.replace(' ', '_')}_${vid}"
    """
    python - ${prefix} "${meta.author_cell_type}" ${h5ad} <<'PY'
import ast
import sys
import anndata
import pandas as pd

prefix, header, h5ad = sys.argv[1], sys.argv[2], sys.argv[3]
log = open('renee_vs_ours_' + prefix + '.log', 'w')
def say(text):
    print(text)
    log.write(text + '\\n')

b = anndata.read_h5ad(h5ad, backed='r')
if b.raw is None:
    sys.exit('STOP: the h5ad has no raw gene table, so no Ensembl ids')
sym = [str(x) for x in b.var_names]
ens = [str(x) for x in b.raw.var_names]
if len(sym) != len(ens):
    sys.exit('STOP: %d gene names and %d Ensembl ids' % (len(sym), len(ens)))
if len(set(sym)) != len(sym):
    sys.exit('STOP: the gene names are not unique')
s2e = dict(zip(sym, ens))
labels = set(b.obs[header].astype(str))
say('h5ad: %d cells, %d genes, %d clusters in %s' % (b.n_obs, b.n_vars, len(labels), header))

# the table of the paper
r = pd.read_csv('in_renee.csv')
for col in ('clusterName', 'f_score', 'NSForest_markers'):
    if col not in r.columns:
        sys.exit('STOP: the results table has no column ' + col)
r['NSForest_markers'] = r['NSForest_markers'].apply(ast.literal_eval)
if r['clusterName'].duplicated().any() or set(r['clusterName']) != labels:
    sys.exit('STOP: the clusters of the results table are not the clusters of the h5ad')
unknown = sorted({g for gs in r['NSForest_markers'] for g in gs if g not in s2e})
if unknown:
    sys.exit('STOP: markers that are not genes of the h5ad: ' + str(unknown[:20]))
say('results table: %d clusters, %d markers' % (len(r), int(r['marker_count'].sum())))

# the check: our combinations result must be the same, cluster by cluster
c = pd.read_csv('in_combinations.csv').set_index('clusterName')
bad = []
for _, row in r.iterrows():
    name = row['clusterName']
    ours = c.loc[name]
    om = ours['markers']
    om = ast.literal_eval(om) if str(om).startswith('[') else str(om).split(', ')
    if list(om) != list(row['NSForest_markers']):
        bad.append('%s: markers %s versus %s' % (name, list(row['NSForest_markers']), list(om)))
    for k in ('f_score', 'recall', 'PPV', 'onTarget'):
        if k in c.columns and k in r.columns and abs(float(row[k]) - float(ours[k])) > 1e-8:
            bad.append('%s: %s %s versus %s' % (name, k, row[k], ours[k]))
    for k in ('TN', 'FP', 'FN', 'TP', 'marker_count'):
        if k in c.columns and k in r.columns and int(row[k]) != int(ours[k]):
            bad.append('%s: %s %s versus %s' % (name, k, row[k], ours[k]))
if bad:
    say('DIFFERENT from our combinations result in %d places, for example:' % len(bad))
    for line in bad[:10]:
        say('  ' + line)
    log.close()
    sys.exit('STOP: the results table differs from our recompute')
say('CHECK: the results table equals our combinations result: same markers, same order, scores within 1e-8, same counts')

# results, markers: the symbols as given, the ENSG ids from the h5ad
res = r.copy()
res.insert(0, 'software_version', '4.0.0')
res.to_csv('results_symbols_' + prefix + '.csv', index=False)
res.to_pickle('results_symbols_' + prefix + '.pkl')
res[['clusterName', 'NSForest_markers', 'f_score']].to_csv('markers_symbols_' + prefix + '.csv', index=False)
res_e = res.copy()
res_e['NSForest_markers'] = res_e['NSForest_markers'].apply(lambda gs: [s2e[g] for g in gs])
res_e.to_csv('results_ensg_' + prefix + '.csv', index=False)
res_e.to_pickle('results_ensg_' + prefix + '.pkl')
res_e[['clusterName', 'NSForest_markers', 'f_score']].to_csv('markers_ensg_' + prefix + '.csv', index=False)
r.rename(columns={'PPV': 'precision'}).to_csv('plots_input_results_' + prefix + '.csv', index=False)

# binary scores and medians: gene x cluster, from the prep step (its row names are the gene symbols)
for kind, src in (('binary_scores', 'in_binary.csv'), ('medians', 'in_medians.csv')):
    df = pd.read_csv(src, index_col=0)
    df.index.name = 'gene'
    missing = [g for g in df.index if g not in s2e]
    if missing:
        sys.exit('STOP: %s: genes that are not genes of the h5ad: %s' % (kind, missing[:20]))
    if set(df.columns) != labels:
        sys.exit('STOP: %s: the columns are not the clusters of %s' % (kind, header))
    df.to_csv(kind + '_symbols_' + prefix + '.csv')
    df.to_pickle(kind + '_symbols_' + prefix + '.pkl')
    e = df.copy()
    e.index = [s2e[g] for g in df.index]
    e.index.name = 'gene'
    e.to_csv(kind + '_ensg_' + prefix + '.csv')
    e.to_pickle(kind + '_ensg_' + prefix + '.pkl')
    say('%s: %d genes x %d clusters, ENSG and symbols' % (kind, df.shape[0], df.shape[1]))
log.close()
PY
    """
}

// dendOrder_author_cell_type.csv of the paper (the order of Figure 5B; column CrossArea_cluster_order) is passed through.
// A label that is not a cluster name only by its letter case (PAX6_1 for Pax6_1) is changed to the cluster name, and
// every change is written in dendOrder_changes_<prefix>.log; any other difference stops the run.
process passthrough_dendorder_process {
    tag "dendorder_${meta.organ}_${meta.first_author}_${meta.year}_${meta.dataset_version_id}"
    container params.nsforest_image
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(dend, stageAs: 'in_dend.csv'), path(renee, stageAs: 'in_renee.csv')

    output:
    tuple val(meta), path("dendOrder_*.csv"), emit: file
    path "dendOrder_changes_*.log",           emit: changes

    script:
    def vid    = meta.dataset_version_id.toString()[-6..-1]
    def prefix = "${meta.organ}_${meta.first_author}_${meta.journal.replace(' ', '_')}_${meta.year}_${meta.author_cell_type.replace(' ', '_')}_${meta.embedding.replace(' ', '_')}_${vid}"
    """
    python - ${prefix} <<'PY'
import sys
import pandas as pd

prefix = sys.argv[1]
col = 'CrossArea_cluster_order'
d = pd.read_csv('in_dend.csv')
names = list(pd.read_csv('in_renee.csv')['clusterName'])
if col not in d.columns:
    sys.exit('STOP: the dendOrder file has no column ' + col)
by_lower = {}
for n in names:
    by_lower.setdefault(n.lower(), []).append(n)
out, changes = [], []
for x in d[col]:
    if x in names:
        out.append(x)
        continue
    m = by_lower.get(x.lower(), [])
    if len(m) != 1:
        sys.exit('STOP: %s is not a cluster name, not even by letter case' % x)
    changes.append((x, m[0]))
    out.append(m[0])
if sorted(out) != sorted(names):
    sys.exit('STOP: the dendOrder file does not list every cluster once')
d[col] = out
d.to_csv('dendOrder_' + prefix + '.csv', index=False)
with open('dendOrder_changes_' + prefix + '.log', 'w') as f:
    f.write('%d of %d labels in %s changed (letter case only):\\n' % (len(changes), len(out), col))
    for a, b in changes:
        f.write('%s -> %s\\n' % (a, b))
print(open('dendOrder_changes_' + prefix + '.log').read())
PY
    """
}

workflow {
    ['datasets_csv', 'organ', 'custom_h5ad', 'silhouette_datasets_csv', 'renee_results', 'dend_order'].each { p ->
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
        .splitCsv(header: true, sep: ',', quote: '"')

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

    // The cluster header is the one published: one copy of every file, none for the subclass header.
    // (The subclass runs above are still needed: their markers are the parents in the combinations.)
    cluster_meta_ch = meta_ch.filter { m -> m.author_cell_type == params.publish_header }
    cluster_h5ad_ch = h5ad_ch.filter { m, f -> m.author_cell_type == params.publish_header }

    // Step 7: the NS-Forest results of the paper, passed through, with the ENSG and the gene-symbol versions of
    // the results, binary scores and medians; our combinations result is the check that they are the same
    renee_in_ch = prep_output_ch.binary_csv
        .filter { m, f -> m.author_cell_type == params.publish_header }
        .join(prep_output_ch.medians_csv)
        .join(h5ad_ch)
        .join(combinations_process.out.results)
        .map { m, binary, medians, h5ad, combos -> tuple(m, file(params.renee_results), binary, medians, h5ad, combos) }
    passthrough_renee_process(renee_in_ch)
    passthrough_dendorder_process(cluster_meta_ch.map { m -> tuple(m, file(params.dend_order), file(params.renee_results)) })
    results_ch = cluster_meta_ch.map { m -> tuple(m, file(params.renee_results)) }

    // plots, from the results of the paper (with PPV called precision)
    if (!params.skip_plots) {
        plots_process(cluster_h5ad_ch.join(passthrough_renee_process.out.plots_results))
    }

    // Step 8: silhouette on EACH supercluster dataset (a cell is compared only with the cells of its own dataset),
    // joined into one set of outputs in the usual format
    sil_rows_ch = Channel
        .fromPath(params.silhouette_datasets_csv)
        .splitCsv(header: true, sep: ',', quote: '"')
    sil_h5ad_ch = cluster_meta_ch.combine(sil_rows_ch).map { m, row ->
        def m5 = m + [dataset_version_id: row.dataset_version_id, supercluster: row.supercluster,
                      order: row.order as int, expected_cells: row.expected_cells, combined: m]
        tuple(m5, file(row.h5ad_url))
    }
    silhouette_output_ch = silhouette_scatter_over_datasets(sil_h5ad_ch)

    viz_summary_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(results_ch)
    )
    viz_distribution_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
    )
    viz_2D_projection_process(cluster_h5ad_ch)

    compute_summary_stats_process(
        silhouette_output_ch.scores
            .join(silhouette_output_ch.cluster_summary)
            .join(silhouette_output_ch.annotation)
            .join(results_ch)
    )

    // Step 9: everything the run writes, for the S3 manifest and for GitHub
    def s3_results_base = workflow.workDir.parent.toUriString() + '/results'

    outputs_ch = Channel.empty()
        .mix(
            dendrogram_output_ch.svg,
            dendrogram_output_ch.cluster_order,
            dendrogram_output_ch.cluster_sizes,
            dendrogram_output_ch.summary,
            cluster_stats_process.out.results,
            cluster_cid_mapping_process.out.results,
            plot_histograms_process.out.histograms,
            passthrough_renee_process.out.files,
            passthrough_dendorder_process.out.file,
            (params.skip_plots ? Channel.empty() : plots_process.out.plots),
            silhouette_output_ch.scores,
            silhouette_output_ch.cluster_summary,
            silhouette_output_ch.annotation,
            viz_summary_process.out.plots,
            viz_distribution_process.out.plots,
            viz_2D_projection_process.out.plots,
            compute_summary_stats_process.out.summary
        )
        .filter { meta, files -> meta.author_cell_type == params.publish_header }

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
