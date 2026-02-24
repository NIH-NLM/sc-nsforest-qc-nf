/**
 * Compute Silhouette Scores
 *
 * Computes silhouette scores for single-cell clusters using embeddings.
 * Applies the same three-stage ontology filter as filter_adata (tissue/disease/age)
 * so that scores are computed only on the filtered normal-adult cell population.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:        Map with organ, first_author, year, author_cell_type, embedding,
 *                  disease, filter
 *   - h5ad:        Path to adata_filtered.h5ad (output of filter_adata_process)
 *   - uberon_json: Path to uberon_{organ}.json — same file used by filter_adata_process.
 *                  Required when meta.filter == "True".
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
 *   - {cluster_header}_silhouette_scores.csv    Per-cell silhouette scores
 *   - {cluster_header}_cluster_summary.csv      Per-cluster mean/median/std
 *   - {cluster_header}_annotation.json          Dataset metadata + filter params recorded
 */
process compute_silhouette_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    shell '/bin/sh'
    
    input:
    tuple val(meta), path(h5ad), path(uberon_json)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_silhouette_scores.csv"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_cluster_summary.csv"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_annotation.json"),
          emit: results
    
    script:
    def filter_flag     = meta.filter == "True" ? "--filter-normal" : "--no-filter-normal"
    def min_age_val     = params.min_age ?: 15
    def min_cluster_val = params.min_cluster_size ?: 5
    """
    docker run ghcr.io/nih-nlm/scsilhouette:1.0 \
        compute-silhouette \
        --h5ad-path ${h5ad} \
        --cluster-header ${meta.author_cell_type} \
        --embedding-key ${meta.embedding} \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        --disease ${meta.disease} \
        ${filter_flag} \
        --uberon ${uberon_json} \
        --min-age ${min_age_val} \
        --min-cluster-size ${min_cluster_val}
    """
}
