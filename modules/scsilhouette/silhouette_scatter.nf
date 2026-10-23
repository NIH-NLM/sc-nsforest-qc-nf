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
