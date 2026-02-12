/**
 * Compute Summary Statistics
 *
 * Creates dataset-level summary statistics from cluster summaries.
 * Computes median-of-medians and other aggregate metrics across all clusters.
 */
process compute_summary_stats_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    container 'ghcr.io/nih-nlm/scsilhouette:1.0'
    publishDir "${params.outdir}", 
        mode: params.publish_mode,
        pattern: "outputs_*/**"
    
    shell '/bin/sh', '-euo', 'pipefail'
    
    input:
    tuple val(meta), 
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)
    
    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_dataset_summary.csv"),
          emit: summary
    
    script:
    """
    cat > summary_stats.py << 'EOF'
import pandas as pd
import numpy as np
import os

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
    'n_clusters': [len(cluster_summary)],
    'n_cells': [cluster_summary['cluster_size'].sum()],
    'median_silhouette': [median_of_medians],
    'mean_silhouette': [mean_of_medians],
    'std_silhouette': [std_of_medians],
    'high_quality_clusters': [high_quality],
    'medium_quality_clusters': [medium_quality],
    'low_quality_clusters': [low_quality],
    'median_fscore': [median_fscore],
    'mean_fscore': [mean_fscore]
})

# Save
os.makedirs("outputs_${meta.organ}_${meta.first_author}_${meta.year}", exist_ok=True)
summary.to_csv("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_dataset_summary.csv", 
               index=False)

print(f"Dataset Summary:")
print(f"  Median of median silhouette scores: {median_of_medians:.3f}")
print(f"  Clusters: {len(cluster_summary)} total")
print(f"    High quality (>= 0.5): {high_quality}")
print(f"    Medium quality (0.25-0.5): {medium_quality}")
print(f"    Low quality (< 0.25): {low_quality}")
if median_fscore:
    print(f"  Median F-score: {median_fscore:.3f}")
EOF

    python3 summary_stats.py
    """
}
