# nsforest_cli/build_symbol_map.py

from pathlib import Path
import pandas as pd
import subprocess
import re

def build_symbol_map_run(
    *,
    gencode_release: int,
    out_csv: Path,
):
    """
    Download and extract ENSG → symbol map from Gencode GTF for given release.
    Produces a CSV with columns: ensg, symbol
    """
    base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    gtf_gz_release = f"release_{gencode_release}/gencode.v{gencode_release}.annotation.gtf.gz"
    gtf_gz_filename = f"gencode.v{gencode_release}.annotation.gtf.gz"
    gtf_filename = gtf_gz_filename.rstrip(".gz")

    if not Path(gtf_gz_filename).exists():
        curl_url = f"{base_url}/{gtf_gz_release}"
        print(f"[info] Downloading {curl_url}")
        subprocess.run(["curl", "-o", gtf_gz_filename, curl_url], check=True)

    if not Path(gtf_filename).exists():
        print(f"[info] Unzipping {gtf_gz_filename}")
        subprocess.run(["gunzip", "-kf", gtf_gz_filename], check=True)

    ensg_to_symbol = {}

    with open(gtf_filename, "r") as f:
        for line in f:
            if line.startswith("#") or "\tgene\t" not in line:
                continue

            fields = line.strip().split("\t")
            attr = fields[-1]

            match_id = re.search(r'gene_id "([^"]+)"', attr)
            match_name = re.search(r'gene_name "([^"]+)"', attr)

            if match_id and match_name:
                ensg = match_id.group(1).split(".")[0]  # strip ENSG version
                symbol = match_name.group(1)
                ensg_to_symbol[ensg] = symbol

    df = pd.DataFrame(list(ensg_to_symbol.items()), columns=["ensg", "symbol"])
    df = df.drop_duplicates(subset="ensg")
    df.to_csv(out_csv, index=False)

    print(f"[done] Wrote {len(df)} unique ENSG→symbol entries to {out_csv}")

    return None
