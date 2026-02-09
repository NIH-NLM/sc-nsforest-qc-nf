/**
 * NSForest Dendrogram Generation Module
 *
 * Generates hierarchical clustering dendrogram showing relationships between
 * cell type clusters based on median gene expression profiles.
 *
 * Corresponds to DEMO_NS-forest_workflow.ipynb: Section 1
 *
 * @input tuple: Dataset metadata and h5ad file
 *   - author_cell_type: Cluster column name in AnnData.obs
 *   - embedding: Embedding key (e.g., X_umap, X_tsne)
 *   - first_author: First author surname
 *   - year: Publication year
 *   - h5ad: Path to h5ad file
 *   - disease: Disease state
 *   - filter: Filter normal cells flag
 *   - tissue: Tissue type
 * @input val organ: Organ/tissue from parameters
 *
 * @output tuple: Results including all output files
 *   - organ: Organ/tissue
 *   - first_author: First author surname
 *   - year: Publication year
 *   - author_cell_type: Cluster column name
 *   - h5ad: Original h5ad file (passed through)
 *   - files: All generated CSV, HTML, and SVG files
 *
 * @output files:
 *   - {cluster_header}_cluster_medians_for_dendrogram.csv: Median expression per cluster
 *   - {cluster_header}_linkage_matrix.csv: Hierarchical clustering linkage matrix
 *   - {cluster_header}_dendrogram.html: Interactive dendrogram visualization
 *   - {cluster_header}_dendrogram.svg: Static dendrogram for publication
 */
process dendrogram_process {
    tag "${organ}_${first_author}_${year}"
    label 'nsforest'
    publishDir "${params.outdir}/${organ}_${first_author}_${year}", mode: params.publish_mode
    
    input:
    tuple val(author_cell_type), val(embedding), val(first_author), val(year), 
          path(h5ad), val(disease), val(filter), val(tissue)
    val organ
    
    output:
    tuple val(organ), val(first_author), val(year), val(author_cell_type), path(h5ad),
          path("outputs_${organ}_${first_author}_${year}/${author_cell_type}_*.{csv,html,svg}")
    
    script:
    """
    nsforest-cli dendrogram \
        --h5ad-path ${h5ad} \
        --cluster-header ${author_cell_type} \
        --organ ${organ} \
        --first-author ${first_author} \
        --year ${year}
    """
}
