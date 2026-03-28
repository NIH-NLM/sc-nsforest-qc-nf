"""
Plot histograms of non-zero median and binary score values.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 histograms

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_hist_nonzero_medians.svg
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_hist_nonzero_binary_scores.svg
"""

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt

from .common_utils import (
    get_output_prefix,
    log_section,
    logger
)


def run_plot_histograms(medians_csv, binary_scores_csv, cluster_header, organ, first_author, journal, year,
                        embedding, dataset_version_id):
    """
    Create histograms of non-zero median and binary score values.
    """
    log_section("NSForest: Plot Histograms")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    logger.info(f"Loading medians: {medians_csv}")
    df_medians = pd.read_csv(medians_csv, index_col=0)

    logger.info(f"Loading binary scores: {binary_scores_csv}")
    df_binary_scores = pd.read_csv(binary_scores_csv, index_col=0)

    # Histogram of non-zero medians
    non_zero_medians = df_medians[df_medians != 0].stack().values
    plt.hist(non_zero_medians, bins=100)
    plt.title("Non-zero medians")
    plt.savefig(f"{prefix}_hist_nonzero_medians.svg")
    plt.close()
    logger.info(f"Saved: {prefix}_hist_nonzero_medians.svg")

    # Histogram of non-zero binary scores
    non_zero_binary_scores = df_binary_scores[df_binary_scores != 0].stack().values
    plt.hist(non_zero_binary_scores, bins=100)
    plt.title("Non-zero binary scores")
    plt.savefig(f"{prefix}_hist_nonzero_binary_scores.svg")
    plt.close()
    logger.info(f"Saved: {prefix}_hist_nonzero_binary_scores.svg")

    logger.info("Histogram plotting complete!")
