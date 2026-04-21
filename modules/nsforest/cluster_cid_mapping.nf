/**
 * Cluster CID Mapping Module
 *
 * Emits a 4-column CSV (cluster_name, skos, manual_mapped_cid, cell_ontology_id)
 * for manual downstream curation of cluster -> cell ontology term assignments.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:  Map with organ, first_author, journal, year, author_cell_type,
 *            embedding, dataset_version_id
 *   - h5ad:  Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit results: tuple(meta, cluster_cid_mapping CSV)
 *   Filename: {organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}_cluster_cid_mapping.csv
 */
process cluster_cid_mapping_process {
    tag "cluster_cid_mapping_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("*.csv"),
          emit: results

    script:
    """
    nsforest-cli cluster-cid-mapping \
        --h5ad-path "${h5ad}" \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}"
    """
}
