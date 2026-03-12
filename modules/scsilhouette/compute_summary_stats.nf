/**
 * Compute Summary Statistics
 *
 * Creates dataset-level summary statistics from cluster summaries.
 * Computes median-of-medians and other aggregate metrics across all clusters.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:             Map with organ, first_author, year, author_cell_type, embedding, doi, etc.
 *   - silhouette_scores:{prefix}_silhouette_scores.csv
 *   - cluster_summary:  {prefix}_cluster_summary.csv
 *   - annotation:       {prefix}_annotation.json
 *   - nsforest_results: {prefix}_results.csv (or NO_FILE sentinel)
 *
 * Output:
 * -------
 * @emit summary: tuple(meta, {prefix}_dataset_summary.csv)
 *   Contains: organ, first_author, year, cluster_header, n_clusters, n_cells,
 *             median/mean/std silhouette, quality tier counts, median/mean F-score,
 *             doi, collection_name, dataset_title, journal
 */
process compute_summary_stats_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta),
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_dataset_summary.csv"),
          emit: summary

    script:
    """
    cat > summary_stats.py << 'EOF'
import pandas as pd
import numpy as np

# Read cluster summary
cluster_summary = pd.read_csv("${cluster_summary}")

# Compute median of medians (primary QC metric)
median_of_medians = cluster_summary['median_silhouette'].median()
mean_of_medians = cluster_summary['median_silhouette'].mean()
std_of_medians = cluster_summary['median_silhouette'].std()

# Count clusters by quality
high_quality = (cluster_summary['median_silhouette'] >= 0.5).sum()
medium_quality = ((cluster_summary['median_silhouette'] >= 0.25) &
                  (cluster_summary['median_silhouette'] < 0.5)).sum()
low_quality = (cluster_summary['median_silhouette'] < 0.25).sum()

# Read NSForest results if available
try:
    nsforest_df = pd.read_csv("${nsforest_results}")
    median_fscore = nsforest_df['f_score'].median()
    mean_fscore = nsforest_df['f_score'].mean()
except:
    median_fscore = None
    mean_fscore = None

# Create summary dataframe
summary = pd.DataFrame({
    'dataset': ["${meta.organ}_${meta.first_author}_${meta.year}"],
    'organ': ["${meta.organ}"],
    'first_author': ["${meta.first_author}"],
    'year': ["${meta.year}"],
    'cluster_header': ["${meta.author_cell_type}"],
    'embedding': ["${meta.embedding}"],
    'doi': ["${meta.doi ?: ''}"],
    'collection_name': ["${meta.collection_name ?: ''}"],
    'dataset_title': ["${meta.dataset_title ?: ''}"],
    'journal': ["${meta.journal ?: ''}"],
    'collection_url': ["${meta.collection_url ?: ''}"],
    'explorer_url': ["${meta.explorer_url ?: ''}"],
    'h5ad_url': ["${meta.h5ad_url ?: ''}"],
    'n_clusters': [len(cluster_summary)],
    'n_cells': [int(cluster_summary['count'].sum())],
    'median_silhouette': [median_of_medians],
    'mean_silhouette': [mean_of_medians],
    'std_silhouette': [std_of_medians],
    'high_quality_clusters': [int(high_quality)],
    'medium_quality_clusters': [int(medium_quality)],
    'low_quality_clusters': [int(low_quality)],
    'median_fscore': [median_fscore],
    'mean_fscore': [mean_fscore]
})

# Save
cluster_header_safe = "${meta.author_cell_type}".replace(" ", "_")
prefix = "${meta.organ}_${meta.first_author}_${meta.year}_" + cluster_header_safe
summary.to_csv(f"{prefix}_dataset_summary.csv", index=False)

print(f"Dataset Summary: ${meta.organ}_${meta.first_author}_${meta.year}")
print(f"  Clusters: {len(cluster_summary)} total")
print(f"  Cells: {int(cluster_summary['count'].sum())}")
print(f"  Median of median silhouette scores: {median_of_medians:.3f}")
print(f"    High quality (>= 0.5): {int(high_quality)}")
print(f"    Medium quality (0.25-0.5): {int(medium_quality)}")
print(f"    Low quality (< 0.25): {int(low_quality)}")
if median_fscore is not None:
    print(f"  Median F-score: {median_fscore:.3f}")
    print(f"  Mean F-score: {mean_fscore:.3f}")
EOF

    python3 summary_stats.py
    """
}
