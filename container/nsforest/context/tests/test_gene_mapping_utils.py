"""
Tests for the gene names, nsforest_cli.gene_mapping_utils.

Each test states its input and what counts as a pass. The tables are made here, small and made up,
so nothing is downloaded. Run inside the image:
    python -m pytest container/nsforest/context/tests
"""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from nsforest_cli import gene_mapping_utils as gm

ROWS = [
    ("ENSG00000000001", "GENEA", "protein_coding", "HGNC:1", True, 48),
    ("ENSG00000000002", "ENSG00000000002", "lncRNA", "", False, 48),   # GENCODE has no symbol: the id is the name
    ("ENSG00000000003", "GENEC-AS1", "lncRNA", "", True, 48),
]


def write_table(folder, release=48, rows=ROWS):
    path = folder / f"gene_mapping_gencode_v{release}.csv"
    pd.DataFrame(rows, columns=["ensembl_gene_id", "gene_name", "gene_type", "hgnc_id",
                                "has_symbol", "gencode_release"]).to_csv(path, index=False)
    return path


def test_every_gene_of_the_table_has_a_name_including_those_without_a_symbol(tmp_path):
    """Input: a table with a gene that has no symbol. Pass: every gene is in the dict, and that gene's name is its id."""
    write_table(tmp_path)
    names = gm.load_gene_mapping(str(tmp_path))
    assert names == {"ENSG00000000001": "GENEA", "ENSG00000000002": "ENSG00000000002",
                     "ENSG00000000003": "GENEC-AS1"}
    assert all(isinstance(name, str) and name for name in names.values())


def test_a_gene_that_is_not_in_the_table_keeps_its_id_as_its_name():
    """Input: markers with an id the table does not have. Pass: the id is the name; no blank, no dropped gene."""
    results = pd.DataFrame({"clusterName": ["c1"], "NSForest_markers": [["ENSG00000000001", "ENSG00000999999"]]})
    _, markers = gm.map_markers_to_symbols(results, {"ENSG00000000001": "GENEA"})
    assert markers["c1"] == ["GENEA", "ENSG00000999999"]


def test_gene_symbols_are_added_to_var_for_every_gene():
    """Input: an AnnData with a known, a symbol-less and an unknown gene. Pass: var['gene_symbol'] has a name for each."""
    adata = anndata.AnnData(X=np.zeros((2, 3), dtype="float32"),
                            var=pd.DataFrame(index=["ENSG00000000001", "ENSG00000000002", "ENSG00000999999"]))
    adata = gm.add_gene_symbols_to_adata(adata, {"ENSG00000000001": "GENEA", "ENSG00000000002": "ENSG00000000002"})
    assert list(adata.var["gene_symbol"]) == ["GENEA", "ENSG00000000002", "ENSG00000999999"]


def test_the_highest_release_in_the_image_is_used(tmp_path):
    """Input: tables of v44 and v48. Pass: v48."""
    write_table(tmp_path, 44)
    write_table(tmp_path, 48)
    assert gm.mapping_file(str(tmp_path)).endswith("gene_mapping_gencode_v48.csv")
    assert gm.release_of("/app/gene_mapping_gencode_v48.csv") == 48


def test_no_table_stops_the_run_instead_of_writing_files_without_names(tmp_path):
    """Input: a folder without a table. Pass: an error that names the problem, not an empty dict."""
    with pytest.raises(FileNotFoundError, match="no gene name table"):
        gm.load_gene_mapping(str(tmp_path))


def test_the_description_for_the_etl_json_names_the_table_its_release_and_checksum(tmp_path):
    """Input: a table and its provenance file. Pass: file, release, rows, SHA-256 and the provenance are recorded."""
    path = write_table(tmp_path)
    (tmp_path / "gene_mapping_gencode_v48.provenance.json").write_text(json.dumps({"source_md5": "abc"}))
    info = gm.describe_mapping(str(tmp_path))
    assert info["file"] == "gene_mapping_gencode_v48.csv" and info["gencode_release"] == 48 and info["rows"] == 3
    assert len(info["sha256"]) == 64 and info["provenance"] == {"source_md5": "abc"}
    assert path.exists()


def test_a_table_without_a_provenance_file_is_described_with_none(tmp_path):
    write_table(tmp_path)
    assert gm.describe_mapping(str(tmp_path))["provenance"] is None
