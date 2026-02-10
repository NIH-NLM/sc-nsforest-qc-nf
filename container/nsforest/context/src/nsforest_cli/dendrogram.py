"""
Generate hierarchical clustering dendrogram.

Corresponds to DEMO_NS-forest_workflow.ipynb: Section 1 - Dendrogram generation
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

from .common_utils import (
    create_output_dir,
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)


def compute_cluster_medians(adata, cluster_header):
    """Compute median expression for each cluster."""
    logger.info("Computing cluster medians for dendrogram...")
    
    clusters = adata.obs[cluster_header].unique()
    n_clusters = len(clusters)
    logger.info(f"Computing medians for {n_clusters} clusters")
    
    median_dict = {}
    for cluster in clusters:
        cluster_mask = adata.obs[cluster_header] == cluster
        cluster_cells = adata[cluster_mask]
        
        if hasattr(adata, 'X') and adata.X is not None:
            median_expr = np.median(
                cluster_cells.X.toarray() if hasattr(cluster_cells.X, 'toarray') else cluster_cells.X,
                axis=0
            )
        else:
            logger.warning("No X matrix found, using raw counts")
            median_expr = np.median(
                cluster_cells.raw.X.toarray() if hasattr(cluster_cells.raw.X, 'toarray') else cluster_cells.raw.X,
                axis=0
            )
        
        median_dict[cluster] = median_expr.flatten()
    
    median_df = pd.DataFrame(median_dict, index=adata.var_names).T
    logger.info(f"Median matrix shape: {median_df.shape}")
    
    return median_df


def generate_dendrogram(median_df, method='ward', metric='euclidean'):
    """Generate hierarchical clustering dendrogram."""
    logger.info(f"Computing hierarchical clustering (method={method}, metric={metric})")
    
    distances = pdist(median_df.values, metric=metric)
    linkage_matrix = linkage(distances, method=method)
    dendro = dendrogram(linkage_matrix, labels=median_df.index.tolist(), no_plot=True)
    
    return linkage_matrix, dendro


def plot_dendrogram_plotly(dendro, cluster_labels, output_prefix):
    """Create interactive dendrogram using Plotly."""
    logger.info("Creating interactive dendrogram plot...")
    
    # Convert to numpy arrays (dendro returns lists)
    icoord = np.array(dendro['icoord'])
    dcoord = np.array(dendro['dcoord'])
    
    fig = go.Figure()
    
    # Add lines for dendrogram
    for i in range(len(icoord)):
        fig.add_trace(go.Scatter(
            x=icoord[i],
            y=dcoord[i],
            mode='lines',
            line=dict(color='black', width=1),
            hoverinfo='skip',
            showlegend=False
        ))
    
    # Generate tick positions (center of each cluster group)
    tick_positions = np.arange(5.0, len(dendro['ivl']) * 10 + 5, 10)
    
    fig.update_layout(
        title="Hierarchical Clustering Dendrogram",
        xaxis=dict(
            title="Cluster",
            tickmode='array',
            tickvals=tick_positions,
            ticktext=dendro['ivl'],
            tickangle=-90
        ),
        yaxis=dict(title="Distance"),
        width=1400,
        height=600,
        showlegend=False,
        margin=dict(b=150)
    )
    
    html_path = f"{output_prefix}_dendrogram.html"
    svg_path = f"{output_prefix}_dendrogram.svg"
    
    pio.write_html(fig, html_path)
    fig.write_image(svg_path)
    
    logger.info(f"Saved: {html_path}")
    logger.info(f"Saved: {svg_path}")


def run_dendrogram(h5ad_path, cluster_header, organ, first_author, year,
                   method='ward', metric='euclidean'):
    """Main function to generate dendrogram."""
    log_section("NSForest: Dendrogram Generation")
    
    output_dir = create_output_dir(organ, first_author, year)
    output_prefix = get_output_prefix(output_dir, cluster_header)
    
    adata = load_h5ad(h5ad_path, cluster_header)
    median_df = compute_cluster_medians(adata, cluster_header)
    
    median_path = f"{output_prefix}_cluster_medians_for_dendrogram"
    median_df.to_csv(f"{median_path}.csv")
    logger.info(f"Saved: {median_path}.csv")
    
    linkage_matrix, dendro = generate_dendrogram(median_df, method=method, metric=metric)
    
    linkage_df = pd.DataFrame(
        linkage_matrix,
        columns=['cluster1', 'cluster2', 'distance', 'n_items']
    )
    linkage_path = f"{output_prefix}_linkage_matrix.csv"
    linkage_df.to_csv(linkage_path, index=False)
    logger.info(f"Saved: {linkage_path}")
    
    plot_dendrogram_plotly(dendro, median_df.index.tolist(), output_prefix)
    
    logger.info("Dendrogram generation complete!")
