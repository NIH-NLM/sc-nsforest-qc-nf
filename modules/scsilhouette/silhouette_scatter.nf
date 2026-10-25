/**
 * Silhouette scores spread over many tasks (replaces compute_silhouette_process for large datasets)
 *
 * The exact silhouette of a cell needs its distance to every other cell. The embedding is small (two
 * columns for a UMAP), so every task reads all of it and scores its own range of cells. Nothing is
 * sampled or approximated: the scores are those of sklearn.metrics.silhouette_samples.
 *
 * The scatter and the gather are done here, in Nextflow. scsilhouette only has the silhouette logic:
 *
 *   silhouette_inputs_process : scsilhouette silhouette-inputs   (reads obs and the embedding, not X)
 *   scatter                   : this workflow cuts 0..n_cells into ranges of params.silhouette_chunk_cells
 *   silhouette_cells_process  : scsilhouette silhouette-cells    (one task per range)
 *   gather                    : silhouette_summary_process joins the score files in cell order with
 *                               shell, then runs scsilhouette silhouette-summary
 *
 * Use:  silhouette_scatter_gather(h5ad_ch)   // h5ad_ch: tuple(meta, h5ad)
 * Emits: scores, cluster_summary, annotation   (tuple(meta, file), as compute_silhouette_process)
 */
process silhouette_inputs_process {
    tag "silhouette_inputs_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.author_cell_type}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("silhouette_embedding_*.npy"),
          path("silhouette_labels_*.npy"),
          path("silhouette_cells_*.csv"),
          path("silhouette_info_*.json"),
          emit: files

    script:
    """
    scsilhouette silhouette-inputs \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --embedding-key "${meta.embedding}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --dataset-version-id "${meta.dataset_version_id}"
    """
}

process silhouette_cells_process {
    tag "silhouette_cells_${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type}_${start}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'

    input:
    tuple val(meta), path(embedding), path(labels), val(start), val(end)

    output:
    tuple val(meta), val(start), path("silhouette_cells_${start}_${end}.csv"), emit: scores

    script:
    """
    scsilhouette silhouette-cells \
        --embedding ${embedding} \
        --labels ${labels} \
        --start ${start} \
        --end ${end} \
        --output silhouette_cells_${start}_${end}.csv
    """
}

process silhouette_summary_process {
    tag "silhouette_summary_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.author_cell_type}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(cells), path(info), path(score_files)    // score_files: in cell order

    output:
    tuple val(meta), path("silhouette_scores*.csv"),            emit: scores
    tuple val(meta), path("scsilhouette_cluster_summary*.csv"), emit: cluster_summary
    tuple val(meta), path("annotation*.json"),                  emit: annotation

    script:
    """
    # join the score files: the header of the first, then the rows of every file, in the order given
    first=1
    for f in ${score_files}; do
        if [ \$first -eq 1 ]; then cat \$f > all_scores.csv; first=0; else tail -n +2 \$f >> all_scores.csv; fi
    done

    scsilhouette silhouette-summary \
        --cells ${cells} \
        --info ${info} \
        --scores all_scores.csv \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --dataset-version-id "${meta.dataset_version_id}" \
        --disease "${meta.disease}"
    """
}

workflow silhouette_scatter_gather {
    take:
    h5ad_ch    // tuple(meta, h5ad)

    main:
    inputs = silhouette_inputs_process(h5ad_ch)

    // scatter: the ranges of cells, from the number of cells in the info file
    ranges_ch = inputs.files.flatMap { meta, embedding, labels, cells, info ->
        def n_cells = new groovy.json.JsonSlurper().parseText(info.text).n_cells as int
        def size    = params.silhouette_chunk_cells as int
        (0..<n_cells).step(size).collect { start ->
            tuple(meta, embedding, labels, start, Math.min(start + size, n_cells))
        }
    }
    cells_out = silhouette_cells_process(ranges_ch)

    // gather: all the score files of a dataset, in cell order
    ordered_ch = cells_out.scores
        .map { meta, start, scores -> tuple(meta, [start, scores]) }
        .groupTuple()
        .map { meta, pairs -> tuple(meta, pairs.sort { it[0] as int }.collect { it[1] }) }

    summary_in = inputs.files
        .map { meta, embedding, labels, cells, info -> tuple(meta, cells, info) }
        .join(ordered_ch)
    summary = silhouette_summary_process(summary_in)

    emit:
    scores          = summary.scores
    cluster_summary = summary.cluster_summary
    annotation      = summary.annotation
}

/**
 * Silhouette scores of SEVERAL datasets (the five CellxGene superclusters), joined into ONE set of outputs
 *
 * Each dataset is scored on its own, as in silhouette_scatter_gather: a cell is compared only with the cells of
 * its own dataset. The cell tables and the score files of the datasets are then stacked, in the order given by
 * meta.order, and the usual silhouette_summary_process writes the usual files (no extra column), once.
 *
 *   meta.combined       the meta of the joined result (the outputs are emitted with it)
 *   meta.supercluster   short name of the dataset (IT, MGE, ...)
 *   meta.order          1, 2, ... the order of the datasets in the joined files
 *   meta.expected_cells the number of cells the dataset must have (the run stops if it differs; blank: no check)
 *
 * Use:  silhouette_scatter_over_datasets(h5ad_ch)   // tuple(meta, h5ad), one for each dataset
 * Emits: scores, cluster_summary, annotation   (tuple(meta.combined, file))
 */
process silhouette_gather_process {
    tag "silhouette_gather_${meta.supercluster}"
    container params.nsforest_image
    cpus 1
    memory '8 GB'

    input:
    tuple val(meta), path(score_files)    // in cell order

    output:
    tuple val(meta), path("all_scores_${meta.supercluster}.csv"), emit: scores

    script:
    """
    # the header of the first file, then the rows of every file, in the order given
    first=1
    for f in ${score_files}; do
        if [ \$first -eq 1 ]; then cat \$f > all_scores_${meta.supercluster}.csv; first=0; else tail -n +2 \$f >> all_scores_${meta.supercluster}.csv; fi
    done
    """
}

process silhouette_concat_process {
    tag "silhouette_concat_${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type}"
    container params.nsforest_image
    cpus 2
    memory '32 GB'

    input:
    tuple val(meta), val(parts_json), path(cells_files), path(info_files), path(score_files)

    output:
    tuple val(meta), path("silhouette_cells_all.csv"), path("silhouette_info_all.json"), path("all_scores_all.csv"), emit: joined

    script:
    """
    cat > parts.json << 'PARTS_EOF'
${parts_json}
PARTS_EOF
    python - ${cells_files} ${info_files} ${score_files} <<'PY'
import json
import sys
import pandas as pd

parts = json.load(open('parts.json'))
n = len(parts)
args = sys.argv[1:]
cells_f, info_f, score_f = args[:n], args[n:2 * n], args[2 * n:3 * n]
if len(args) != 3 * n:
    sys.exit('STOP: %d files for %d datasets' % (len(args), n))

cells, infos, scores = [], [], []
for p, cf, jf, sf in zip(parts, cells_f, info_f, score_f):
    c = pd.read_csv(cf, index_col=0, dtype=str)
    info = json.load(open(jf))
    s = pd.read_csv(sf)['silhouette_score']
    if not (len(c) == len(s) == info['n_cells']):
        sys.exit('STOP: %s: %d cells, %d scores, %d in the info file' % (p['supercluster'], len(c), len(s), info['n_cells']))
    exp = str(p.get('expected_cells') or '').strip()
    if exp and int(exp) != info['n_cells']:
        sys.exit('STOP: %s has %d cells, expected %s' % (p['supercluster'], info['n_cells'], exp))
    print('%s: %d cells, %d clusters, scores %.3f to %.3f' % (p['supercluster'], info['n_cells'], info['n_clusters'], s.min(), s.max()))
    cells.append(c)
    infos.append(info)
    scores.append(s)

base = infos[0]
for p, info in zip(parts[1:], infos[1:]):
    for k in ('cluster_header', 'embedding_key', 'embedding_dtype'):
        if info[k] != base[k]:
            sys.exit('STOP: %s: %s is %s, the first dataset has %s' % (p['supercluster'], k, info[k], base[k]))

allc = pd.concat(cells)
if allc.index.duplicated().any():
    sys.exit('STOP: %d cell names occur in more than one dataset' % int(allc.index.duplicated().sum()))
alls = pd.concat(scores, ignore_index=True)
if len(allc) != len(alls):
    sys.exit('STOP: %d cells and %d scores' % (len(allc), len(alls)))

out = dict(base)
out['n_cells'] = int(len(allc))
out['n_clusters'] = int(allc[base['cluster_header']].nunique())
out['obs_keys'] = [k for k in base['obs_keys'] if all(k in i['obs_keys'] for i in infos)]
out['obsm_keys'] = [k for k in base['obsm_keys'] if all(k in i['obsm_keys'] for i in infos)]
allc.to_csv('silhouette_cells_all.csv')
alls.to_frame('silhouette_score').to_csv('all_scores_all.csv', index=False)
json.dump(out, open('silhouette_info_all.json', 'w'), indent=2)
print('joined: %d cells, %d clusters' % (out['n_cells'], out['n_clusters']))
PY
    """
}

workflow silhouette_scatter_over_datasets {
    take:
    h5ad_ch    // tuple(meta, h5ad), one for each dataset

    main:
    inputs = silhouette_inputs_process(h5ad_ch)

    ranges_ch = inputs.files.flatMap { meta, embedding, labels, cells, info ->
        def n_cells = new groovy.json.JsonSlurper().parseText(info.text).n_cells as int
        def size    = params.silhouette_chunk_cells as int
        (0..<n_cells).step(size).collect { start ->
            tuple(meta, embedding, labels, start, Math.min(start + size, n_cells))
        }
    }
    cells_out = silhouette_cells_process(ranges_ch)

    ordered_ch = cells_out.scores
        .map { meta, start, scores -> tuple(meta, [start, scores]) }
        .groupTuple()
        .map { meta, pairs -> tuple(meta, pairs.sort { it[0] as int }.collect { it[1] }) }
    gathered = silhouette_gather_process(ordered_ch)

    // all the datasets of one result, in the order of meta.order
    joined_in = inputs.files
        .map { meta, embedding, labels, cells, info -> tuple(meta, cells, info) }
        .join(gathered)
        .map { meta, cells, info, scores ->
            tuple(meta.combined, meta.order as int, meta.supercluster, meta.expected_cells, cells, info, scores)
        }
        .groupTuple()
        .map { combined, orders, names, expected, cells, infos, scores ->
            def idx   = (0..<orders.size()).toList().sort { orders[it] }
            def parts = idx.collect { [supercluster: names[it], expected_cells: expected[it]] }
            tuple(combined, groovy.json.JsonOutput.toJson(parts),
                  idx.collect { cells[it] }, idx.collect { infos[it] }, idx.collect { scores[it] })
        }
    joined = silhouette_concat_process(joined_in)
    summary = silhouette_summary_process(joined.joined.map { meta, cells, info, scores -> tuple(meta, cells, info, scores) })

    emit:
    scores          = summary.scores
    cluster_summary = summary.cluster_summary
    annotation      = summary.annotation
}
