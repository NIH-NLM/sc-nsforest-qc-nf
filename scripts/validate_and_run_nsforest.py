import anndata as ad
from anndata import AnnData

def filter_adata(
    adata: AnnData,
    obs_key: str,
    values: list[str],
    mode: str = "exact",              # "exact", "contains", or "regex"
    na_policy: str = "drop",          # "drop", "keep", or "match"
    case_insensitive: bool = True,
    invert: bool = False
) -> AnnData:
    import os
    import sys
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import matplotlib.pyplot as plt
    # Insert local NSForest path
    sys.path.insert(0, os.path.abspath("/home/jovyan/session_data/NSForest"))
    import nsforest as ns
    import pandas as pd
    import numpy as np
    from anndata import AnnData
    
    if obs_key not in adata.obs.columns:
        raise KeyError(f"obs_key '{obs_key}' not in adata.obs")

    s = adata.obs[obs_key].astype("string")

    if na_policy == "drop":
        mask_valid = s.notna()
        s = s[mask_valid]
    elif na_policy == "keep":
        mask_valid = pd.Series(True, index=s.index)
    elif na_policy == "match":
        s = s.fillna("NA")
        mask_valid = pd.Series(True, index=s.index)
    else:
        raise ValueError("na_policy must be 'drop' | 'keep' | 'match'")

    if case_insensitive:
        s_cmp = s.str.lower()
        vals = [v.lower() for v in values]
    else:
        s_cmp = s
        vals = values

    if mode == "exact":
        mask_match = pd.Series(False, index=s_cmp.index)
        for v in vals:
            mask_match |= (s_cmp == v)
    elif mode == "contains":
        mask_match = pd.Series(False, index=s_cmp.index)
        for v in vals:
            mask_match |= s_cmp.str.contains(v, regex=False)
    elif mode == "regex":
        pattern = "(" + "|".join(vals) + ")"
        mask_match = s_cmp.str.contains(pattern, regex=True)
    else:
        raise ValueError("mode must be 'exact' | 'contains' | 'regex'")

    final_mask = mask_valid & ( ~mask_match if invert else mask_match )
    final_mask_np = final_mask.to_numpy(dtype=bool)
    return adata[final_mask_np].copy()

def validate_and_run_nsforest(csv_path, output_path=None):
    import os
    import sys
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import matplotlib.pyplot as plt
    # Insert local NSForest path
    sys.path.insert(0, os.path.abspath("/home/jovyan/session_data/NSForest"))
    import nsforest as ns
    import pandas as pd
    import numpy as np
    from anndata import AnnData
    
    df = pd.read_csv(csv_path)
    outputs = []

    for idx, row in df.iterrows():

        # 0. Set up
        organ     = row.get('tissue', None)
        author    = row.get('author', f'Entry {idx}')
        year      = row.get('publication_date', None)
        h5ad      = row.get('h5ad_file', None)
        label_key = row.get('label_key', None)

        output_folder = f"outputs_{organ}_{author}_{year}/"
        os.makedirs(output_folder, exist_ok=True)

        cluster_header = label_key
        outputfilename_suffix = cluster_header
        outputfilename_prefix = cluster_header 
        cluster_header = label_key

        # 1. Data 
        if not h5ad or not label_key:
            outputs.append(f"Skipping {author}: missing h5ad or label_key\n")
            continue

        outputs.append(f"\n=== {author} ===")
        outputs.append(f"Loading: {h5ad}")

        try:
            if h5ad.startswith("http"):
                outputs.append(f"Skipped: Remote .h5ad file detected (not downloaded)")
                continue

            adata = sc.read_h5ad(h5ad)

            # Print and save adata summary
            adata_summary = str(adata)
            print("pre-filtering summary")
            print(adata_summary)
            outputs.append(f"\n{adata_summary}\n")

            if cluster_header not in adata.obs:
                outputs.append(f"Missing cluster header: {cluster_header}")
                continue

            # Filter for only normal samples
            adata_filtered = filter_adata(adata, obs_key="disease", values=["normal"])
            outputs.append(f"Original cells: {adata.shape[0]}")
            outputs.append(f"Filtered normal cells: {adata_filtered.shape[0]}")

            # Then filter by tissue
            adata_filtered_tissue = filter_adata(adata_filtered, obs_key="tissue", values=[str(organ)])
            outputs.append(f"Normal cells: {adata_filtered.shape[0]}")
            outputs.append(f"Filtered tissue cells: {adata_filtered_tissue.shape[0]}")

            adata = adata_filtered_tissue
            
            # Print and save adata summary after filtering
            adata_summary = str(adata)
            print(adata_summary)
            outputs.append(f"\n{adata_summary}\n")

            if cluster_header not in adata.obs:
                outputs.append(f"Missing cluster header: {cluster_header}")
                continue

#            output_folder = os.path.dirname(h5ad)
            filename = os.path.basename(h5ad)
            outputfilename_prefix = filename.replace(".h5ad", "")
            outputfilename_suffix = outputfilename_prefix

            # 2. Clusters
            # number of clusters
            n_clusters = adata.obs[cluster_header].nunique()
            print(f"Number of clusters: {n_clusters}")
            outputs.append(f"Number of clusters: {n_clusters}\n")

            ## auto-adjust figsize
            fig_width = int(n_clusters/5)
            fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

            ## dendrogram and save svg
            ns.pp.dendrogram(
                adata,
                cluster_header,
                figsize = (fig_width, fig_height),
                tl_kwargs = {'optimal_ordering': True},
                save = "svg",
                output_folder = output_folder,
                outputfilename_suffix = outputfilename_suffix
            )

            # cluster sizes
            df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
            print(str(df_cluster_sizes))
            outputs.append(str(df_cluster_sizes))

            # save 
            df_cluster_sizes.to_csv(
                os.path.join(
                    output_folder,
                    outputfilename_prefix + "_cluster_sizes.csv"
                )
            )

             # cluster order
            cluster_order = [
                x.strip() for x in adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
            ]

            # save
            pd.DataFrame({"cluster_order": cluster_order}).to_csv(
                os.path.join(
                    output_folder,
                    outputfilename_prefix + "_cluster_order.csv"),
                index=False
            )

            # summary statistics of data (normal cells)
            df_normal = pd.DataFrame(
                {"n_obs": [adata.n_obs],
                 "n_vars": [adata.n_vars],
                 "n_clusters": [n_clusters]}
            )
            print(str(df_normal))            
            outputs.append(str(df_normal))

            # save
            df_normal.to_csv(
                os.path.join(
                    output_folder, outputfilename_prefix + "_summary_normal.csv"
                ),
                index=False
            )

            # 3. NS-Forest
            # prep
            ## make a copy b/c the median step will only keep the positive genes
            ## keep the original data for plotting
            adata_prep = adata.copy()
            
            ## get medians
            adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)
            
            ## get binary scores
            adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)
            
            ## check medians
            df_medians = adata_prep.varm['medians_' + cluster_header]
            print(df_medians.shape)
            df_medians.head()
            
            ## check binary scores
            df_binary_scores = adata_prep.varm['binary_scores_' + cluster_header]
            print(df_binary_scores.shape)
            df_binary_scores.head()
            
            ## save csv and pkl
            df_medians.to_csv(output_folder + outputfilename_prefix + "_medians.csv")
            df_medians.to_pickle(output_folder + outputfilename_prefix + "_medians.pkl")
            
            df_binary_scores.to_csv(output_folder + outputfilename_prefix + "_binary_scores.csv")
            df_binary_scores.to_pickle(output_folder + outputfilename_prefix + "_binary_scores.pkl")

            # histograms of non-zero values [TO-DO: nice to have them as functions -- Beverly]
            non_zero_medians = df_medians[df_medians != 0].stack().values
            plt.hist(non_zero_medians, bins=100)
            plt.title("Non-zero medians")
            
            plt.savefig(output_folder + "hist_nonzero_medians_" + outputfilename_suffix + ".svg")
            plt.show()

            non_zero_binary_scores = df_binary_scores[df_binary_scores != 0].stack().values
            plt.hist(non_zero_binary_scores, bins=100)
            plt.title("Non-zero binary scores")
            
            plt.savefig(output_folder + "hist_nonzero_binary_scores_" + outputfilename_suffix + ".svg")
            plt.show()

            # run NSForest()
            results = ns.nsforesting.NSForest(
                adata_prep,
                cluster_header,
                save = True,
                save_supplementary = True,
                output_folder = output_folder,
                outputfilename_prefix = outputfilename_prefix
            )

            results

            ## 4. Plotting
            #### **load NS-Forest results** (copy set up and load pkl)

            ## set results to plot
            results_to_plot = results

            #### boxplots
            ns.pl.boxplot(results_to_plot, "f_score", save = "html", output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.boxplot(results_to_plot, "precision", save = "html", output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.boxplot(results_to_plot, "recall", save = "html", output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.boxplot(results_to_plot, "recall", save = "html", output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.boxplot(results_to_plot, "onTarget", save = "html", output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            
            # scatter plots w.r.t. cluster
            ns.pl.scatter_w_clusterSize(results, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.scatter_w_clusterSize(results, "precision", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.scatter_w_clusterSize(results, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
            ns.pl.scatter_w_clusterSize(results, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

            # #### scanpy plots [TO-DO: map ENSGs to gene symbols using Ray's repo -- Anne & Ray]
            ## make the markers dictionary for plotting
            markers_dict = dict(zip(results_to_plot["clusterName"], results_to_plot["NSForest_markers"]))
            
            ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False,
                          save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
            ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False, standard_scale = 'var',
                          save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix + "_scaled")
            ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False,
                                save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
            ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False, standard_scale = 'var',
                                save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix + "_scaled")
            ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False,
                             save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
            ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = True, use_raw = False, standard_scale = 'var',
                             save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix + "_scaled")
            
        except Exception as e:
            outputs.append(f"Error processing {author}: {str(e)}")

    if output_path:
        with open(output_path, "w") as f:
            f.write("\n".join(outputs))
    else:
        print("\n".join(outputs))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Validate h5ad and run nsforest")
    parser.add_argument("csv", help="CSV file with h5ad info")
    parser.add_argument("output", nargs="?", default=None, help="Optional: write output to file")
    args = parser.parse_args()

    validate_and_run_nsforest(args.csv, args.output)
