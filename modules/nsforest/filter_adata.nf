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
 *   - meta:         Map with organ, first_author, journal, year, author_cell_type, embedding, vid, filter
 *   - h5ad:         Path to input h5ad (downloaded from CellxGene)
 *   - uberon_json:  Path to uberon_{organ}.json from cellxgene-harvester resolve-uberon
 *   - disease_json: Path to disease_normal.json from cellxgene-harvester resolve-disease
 *   - hsapdv_json:  Path to hsapdv_adult_N.json from cellxgene-harvester resolve-hsapdv
 *
 * Output:
 * -------
 * @emit results: tuple(meta, adata_filtered.h5ad, [stats CSVs and SVGs])
 *   Flat filenames: {organ}_{first_author}_{journal}_{year}_{cluster_header_safe}_{embedding_safe}_{vid}_adata_filtered.h5ad
 *                   {organ}_{first_author}_{journal}_{year}_{cluster_header_safe}_{embedding_safe}_{vid}*.{csv,svg}
 */
process filter_adata_process {
    tag "filter_adata_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "*.{h5ad,csv,svg,log}"
    input:
    tuple val(meta), path(h5ad)
    path (uberon_json,  stageAs: 'uberon.json')
    path (disease_json, stageAs: 'disease.json')
    path (hsapdv_json,  stageAs: 'hsapdv.json')

    output:
    tuple val(meta), path("*adata_filtered.h5ad"),               emit: h5ad
    tuple val(meta), path("*_cluster_sizes_before_filter.csv"),  emit: cluster_sizes
    tuple val(meta), path("*_cluster_order_before_filter.csv"),  emit: cluster_order
    tuple val(meta), path("*_summary_before_filter.csv"),        emit: summary
    tuple val(meta), path("*_before_filter.svg"),                emit: svg
    

    script:
    def filter_flag     = meta.filter == "True" ? "--filter-normal" : ""
    def obs_col_flag    = meta.filter_obs_column ? "--filter-obs-column ${meta.filter_obs_column}" : ""
    def obs_val_flag    = meta.filter_obs_value  ? "--filter-obs-value ${meta.filter_obs_value}"   : ""
    def min_cluster_val = params.min_cluster_size ?: 5
    """
    nsforest-cli filter-adata \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
        ${filter_flag} \
        --uberon ${uberon_json} \
        --disease ${disease_json} \
        --hsapdv ${hsapdv_json} \
        --min-cluster-size ${min_cluster_val} \
        --embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}" \
	${obs_col_flag} \
        ${obs_val_flag}

    """
}
