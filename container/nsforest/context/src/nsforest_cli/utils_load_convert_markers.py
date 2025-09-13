import pandas as pd
import ast
import anndata as ad

def load_and_convert_markers(results_csv,
                             symbol_map_csv):

    """
    Load NS-Forest marker results and return a markers_dict
    with gene symbols (converted from Ensembl IDs in list form).
    """
    cluster_col='clusterName'
    markers_col='NSForest_markers'
    ensembl_col='ensg'
    symbol_col='symbol'

    
    # Load results
    df = pd.read_csv(results_csv)

    # Load symbol map and create dictionary
    symbol_map_df = pd.read_csv(symbol_map_csv)
    symbol_map = dict(zip(symbol_map_df[ensembl_col], symbol_map_df[symbol_col]))

    # Parse list strings into actual lists
    df[markers_col] = df[markers_col].apply(ast.literal_eval)

    # Map each Ensembl ID in the list to gene symbol
    df[symbol_col] = df[markers_col].apply(
        lambda id_list: [symbol_map.get(gid, gid) for gid in id_list]
    )

    # Build markers_dict: {cluster_name: [gene1, gene2, ...]}
    markers_dict = df.set_index(cluster_col)[symbol_col].to_dict()

    return markers_dict


def convert_adata_varnames_with_symbol_map(
    adata: ad.AnnData,
    symbol_map_csv_path: str,
) -> ad.AnnData:
    """
    Convert adata.var_names from ENSG IDs to gene symbols using a provided symbol map.
    The original ENSG IDs are saved in adata.var['orig_var_names'].

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with var_names as ENSG IDs.
    symbol_map_csv_path : Optional[str]
        Path to CSV containing ENSG ID -> gene symbol mapping.
        Assumes first column is ENSG ID, second column is gene symbol.

    Returns
    -------
    AnnData
        Updated AnnData object with var_names replaced by gene symbols.
    """
    if symbol_map_csv_path is None:
        return adata

    symbol_map_df = pd.read_csv(symbol_map_csv_path)
    symbol_map = dict(zip(symbol_map_df.iloc[:, 0], symbol_map_df.iloc[:, 1]))

    adata.var["orig_var_names"] = adata.var_names
    adata.var_names = [symbol_map.get(gene, gene) for gene in adata.var_names]

    return adata

