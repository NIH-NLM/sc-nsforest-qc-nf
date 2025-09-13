import pandas as pd

def load_symbol_map(symbol_map_csv: str | None) -> dict[str, str]:
    if not symbol_map_csv:
        return {}

    df = pd.read_csv(symbol_map_csv)
    df = df.dropna().drop_duplicates(subset=["ensg"])

    # Normalize ENSG keys (strip version suffixes if needed)
    ensg_to_symbol = {}
    for ensg, symbol in zip(df["ensg"], df["symbol"]):
        if "." in ensg:
            ensg = ensg.split(".")[0]
        ensg_to_symbol[ensg] = symbol

    return ensg_to_symbol

def apply_symbol_map(gene_ids: list[str], symbol_map: dict[str, str]) -> list[str]:
    """Map list of gene IDs to symbols using provided mapping"""
    return [symbol_map.get(g.split(".")[0], g) for g in gene_ids]
