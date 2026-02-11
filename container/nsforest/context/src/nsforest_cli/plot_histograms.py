"""
Plot histograms of non-zero median and binary score values.

Creates two SVG histograms showing the distribution of non-zero values.
"""

import pandas as pd
import matplotlib.pyplot as plt

from .common_utils import (
    create_output_dir,
    log_section,
    logger
)


def run_plot_histograms(medians_csv, binary_scores_csv, cluster_header, organ, first_author, year):
    """
    Create histograms of non-zero values.
    
    Reads:
    - medians.csv
    - binary_scores.csv
    
    Creates:
    - hist_nonzero_medians_{cluster_header}.svg
    - hist_nonzero_binary_scores_{cluster_header}.svg
    """
    log_section("NSForest: Plot Histograms")
    
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_suffix = cluster_header
    
    # Load data
    logger.info(f"Loading medians: {medians_csv}")
    df_medians = pd.read_csv(medians_csv, index_col=0)
    
    logger.info(f"Loading binary scores: {binary_scores_csv}")
    df_binary_scores = pd.read_csv(binary_scores_csv, index_col=0)
    
    # Histogram of non-zero medians
    logger.info("Creating histogram of non-zero medians...")
    non_zero_medians = df_medians[df_medians != 0].stack().values
    
    plt.figure(figsize=(10, 6))
    plt.hist(non_zero_medians, bins=100)
    plt.title("Non-zero medians")
    plt.xlabel("Median expression")
    plt.ylabel("Frequency")
    
    hist_medians_path = output_folder + "/" + "hist_nonzero_medians_" + outputfilename_suffix + ".svg"
    plt.savefig(hist_medians_path)
    plt.close()
    logger.info(f"Saved: hist_nonzero_medians_{outputfilename_suffix}.svg")
    
    # Histogram of non-zero binary scores
    logger.info("Creating histogram of non-zero binary scores...")
    non_zero_binary_scores = df_binary_scores[df_binary_scores != 0].stack().values
    
    plt.figure(figsize=(10, 6))
    plt.hist(non_zero_binary_scores, bins=100)
    plt.title("Non-zero binary scores")
    plt.xlabel("Binary score")
    plt.ylabel("Frequency")
    
    hist_binary_path = output_folder + "/" + "hist_nonzero_binary_scores_" + outputfilename_suffix + ".svg"
    plt.savefig(hist_binary_path)
    plt.close()
    logger.info(f"Saved: hist_nonzero_binary_scores_{outputfilename_suffix}.svg")
    
    logger.info("Histogram plotting complete!")
