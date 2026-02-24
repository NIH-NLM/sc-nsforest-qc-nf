/**
 * Filter AnnData Module
 *
 * Three-stage filtering using ontology IDs — identical logic to cellxgene-harvester:
 *   1. Tissue  — tissue_ontology_term_id.isin(obo_ids) from uberon_{organ}.json
 *   2. Disease — disease_ontology_term_id == PATO:0000461
 *   3. Age     — development_stage text parsed, keeps cells >= min_age years
 * Then:
 *   4. Min cluster size — drops clusters with < 5 cells
 *
 * Creates before/after dendrograms and cluster statistics.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map with organ, first_author, year, author_cell_type, filter
 *   - h5ad:        Path to input h5ad (from CellxGene, must have tissue_ontology_term_id
 *                  and disease_ontology_term_id in obs)
 *   - uberon_json: Path to uberon_{organ}.json produced by cellxgene-harvester resolve-uberon.
 *                  Contains obo_ids covering the full anatomical unit (e.g. heart + pericardium
 *                  + all descendant UBERON terms). Same file used for all h5ad in this run.
 *
 * Key params:
 * -----------
 * params.min_age          Minimum donor age in years (default 15)
 * params.min_cluster_size Minimum cells per cluster (default 5)
 * params.publish_mode     Nextflow publishDir mode
 *
 * Output:
 * -------
 * @return tuple:
 *   - meta
 *   - adata_filtered.h5ad
 *   - before/after cluster stats, sizes, cluster order CSVs and dendrograms (SVGs)
 */
process filter_adata_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    input:
    tuple val(meta), path(h5ad), path(uberon_json)
    
    output:
    tuple val(meta), 
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_filtered.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}*.{csv,svg}"),
          emit: results
    
    script:
    def filter_flag     = meta.filter == "True" ? "--filter-normal" : ""
    def min_age_val     = params.min_age ?: 15
    def min_cluster_val = params.min_cluster_size ?: 5
    """
    nsforest-cli filter-adata \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        ${filter_flag} \
        --uberon ${uberon_json} \
        --min-age ${min_age_val} \
        --min-cluster-size ${min_cluster_val}
    """
}
