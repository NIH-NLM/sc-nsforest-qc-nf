/**
 * Filter AnnData Module
 *
 * Three-stage filtering using ontology IDs — identical logic to cellxgene-harvester:
 *   1. Tissue  — tissue_ontology_term_id.isin(obo_ids) from uberon_{organ}.json
 *   2. Disease — disease_ontology_term_id.isin(obo_ids) from disease_normal.json
 *   3. Age     — development_stage_ontology_term_id.isin(obo_ids) from hsapdv_adult_N.json
 * Then:
 *   4. Min cluster size — drops clusters with < min_cluster_size cells
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:         Map with organ, first_author, year, author_cell_type, filter
 *   - h5ad:         Path to input h5ad (downloaded from CellxGene)
 *   - uberon_json:  Path to uberon_{organ}.json from cellxgene-harvester resolve-uberon
 *   - disease_json: Path to disease_normal.json from cellxgene-harvester resolve-disease
 *   - hsapdv_json:  Path to hsapdv_adult_N.json from cellxgene-harvester resolve-hsapdv
 *
 * Output:
 * -------
 * @return tuple:
 *   - meta
 *   - adata_filtered.h5ad
 *   - before/after cluster stats, sizes, cluster order CSVs and dendrograms (SVGs)
 */
process filter_adata_process {
    tag "\"${meta.organ}\"_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "outputs_*/**"

    input:
    tuple val(meta), path(h5ad), path(uberon_json), path(disease_json), path(hsapdv_json)

    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_filtered.h5ad"),
	  path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/*.{csv,svg}",optional: true),
          emit: results

    script:
    def filter_flag     = meta.filter == "True" ? "--filter-normal" : ""
    def min_cluster_val = params.min_cluster_size ?: 5
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
        --min-cluster-size ${min_cluster_val}
    """
}
