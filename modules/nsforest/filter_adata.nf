/**
 * Filter AnnData Module
 *
 * Three-stage filtering using per-row ontology term IDs from harvester CSV:
 *   1. Tissue  — tissue_ontology_term_id
 *   2. Disease — disease_ontology_term_id
 *   3. Age     — development_stage_ontology_term_id
 *   4. Min cluster size
 */
process filter_adata_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), path(uberon_json), path(disease_json), path(hsapdv_json)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_adata_filtered.h5ad"),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}*.{csv,svg}"),
          emit: results

    script:
    def filter_flag     = meta.filter == "True" ? "--filter-normal" : ""
    def min_cluster_val = params.min_cluster_size ?: 5
    def tissue_ids      = meta.tissue_ontology_term_id             ? "--tissue-ontology-term-id \"${meta.tissue_ontology_term_id}\""                         : ""
    def disease_ids     = meta.disease_ontology_term_id            ? "--disease-ontology-term-id \"${meta.disease_ontology_term_id}\""                        : ""
    def hsapdv_ids      = meta.development_stage_ontology_term_id  ? "--development-stage-ontology-term-id \"${meta.development_stage_ontology_term_id}\""   : ""
    """
    nsforest-cli filter-adata \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        ${filter_flag} \
        --uberon ${uberon_json} \
        --disease ${disease_json} \
        --hsapdv ${hsapdv_json} \
        --min-cluster-size ${min_cluster_val} \
        ${tissue_ids} \
        ${disease_ids} \
        ${hsapdv_ids}
    """
}
