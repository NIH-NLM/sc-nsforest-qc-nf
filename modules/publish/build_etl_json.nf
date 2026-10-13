/**
 * Build the master dataset summary JSON (the one file per dataset that the ETL ingests)
 *
 * Writes the one JSON file for a dataset that goes on to the ETL (master_dataset_summary_*.json). It holds the whole
 * cellxgene-harvester-nf record (nothing is left out), the sc-nsforest-qc-nf run
 * (parameters, the dataset summary, the result files) and where the final filtered
 * h5ad file is. The location is in "s3_filtered_h5ad": with --s3_h5ad_prefix it is
 * <prefix>/<file name>; without it (a test) it is the local file name.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:           Map with organ, first_author, year, embedding, dataset_version_id, ...
 *   - summary_csvs:   the CSVs from compute_summary_stats (master_dataset_summary_*.csv is used)
 *   - harvester_json: <organism>_<organ>_harvester_final.json from cellxgene-harvester-nf (the list of
 *                     records of all datasets; the record of meta.dataset_id is used)
 *   - h5ad:           the final filtered h5ad (adata_filtered_*.h5ad); its CellxGene schema version
 *                     (uns['schema_version']) is recorded in sc_nsforest_qc.gene_annotation
 *   - positive_csv:   binary_positive_genes_*.csv (every gene with a binary score above 0, per cluster)
 *   - top10_csv:      binary_genes_top10_*.csv (the top binary genes per cluster from NSForest, with
 *                     the random forest importance and the cluster median)
 *   - uberon_json, disease_json, hsapdv_json, assay_json: the resolve files given to the workflow;
 *     their names and sha256 are recorded in sc_nsforest_qc.resolve_files
 *
 * Output:
 * -------
 * @emit json: tuple(meta, master_dataset_summary_{prefix}.json)
 */
process build_etl_json_process {
    tag "etl_json_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'publish_json'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(summary_csvs), path(harvester_json), path(h5ad), path(positive_csv), path(top10_csv),
          path(input_info), path(sizes_before), path(summary_before)
    path(uberon_json,  stageAs: 'resolve/uberon.json')
    path(disease_json, stageAs: 'resolve/disease.json')
    path(hsapdv_json,  stageAs: 'resolve/hsapdv.json')
    path(assay_json,   stageAs: 'resolve/assay.json')
    path(harvester_csv, stageAs: 'input/harvester.csv')

    output:
    tuple val(meta), path("master_dataset_summary_*.json"), emit: json

    script:
    def h5ad_name = h5ad.name
    def run = groovy.json.JsonOutput.toJson([
        workflow_version:      workflow.manifest.version ?: '',
        session_id:            meta.session_id,
        min_cluster_size:      params.min_cluster_size,
        n_trees:               params.n_trees,
        max_cells_per_cluster: params.max_cells_per_cluster,
        nsforest_seed:         params.nsforest_seed,
        batch_size:            params.batch_size,
        filter_obs_column:     meta.filter_obs_column,
        filter_obs_value:      meta.filter_obs_value,
        author_cell_type:      meta.author_cell_type,
        embedding:             meta.embedding,
        input_filtered_h5ad_dir: params.h5ad_dir ? params.h5ad_dir.toString() : null,
        curation:              [reference: meta.reference, author_cell_type: meta.author_cell_type, embedding: meta.embedding],
        resolve_file_names:    [uberon: file(params.uberon_json).name, disease: file(params.disease_json).name, hsapdv: file(params.hsapdv_json).name, assay: file(params.assay_json).name],
        filtered_h5ad:         h5ad_name,
        s3_filtered_h5ad:      params.s3_h5ad_prefix ? "${params.s3_h5ad_prefix.toString().replaceAll('/+\$', '')}/${h5ad_name}" : h5ad_name,
    ])
    """
    cat > run.json <<'RUN_JSON'
${run}
RUN_JSON
python3 - <<'PY'
import csv, glob, hashlib, json, os
record = next(r for r in json.load(open("${harvester_json}")) if r["dataset"]["dataset_id"] == "${meta.dataset_id}")
run = json.load(open("run.json"))
names = run.pop("resolve_file_names")
run["resolve_files"] = {
    key: {"file": names[key], "sha256": hashlib.sha256(open("resolve/" + key + ".json", "rb").read()).hexdigest()}
    for key in ("uberon", "disease", "hsapdv", "assay")
}
def gene_annotation():
    # the CellxGene schema version of the h5ad this workflow read, and which gene lookup was used
    import h5py
    with h5py.File("${h5ad}", "r") as f:
        version = f["uns"]["schema_version"][()] if "schema_version" in f["uns"] else None
    version = version.decode() if isinstance(version, bytes) else version
    from nsforest_cli.gene_mapping_utils import describe_mapping
    info = describe_mapping()
    return {"h5ad_schema_version": version, "gene_mapping": info}

summary = {}
for f in glob.glob("master_dataset_summary_*.csv"):
    row = next(csv.DictReader(open(f)))
    for k, v in row.items():
        try:
            summary[k] = int(v)
        except ValueError:
            try:
                summary[k] = float(v)
            except ValueError:
                summary[k] = v
def table(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k in ("binary_score", "rf_feature_importance", "cluster_median"):
            if k in r and r[k] != "":
                r[k] = float(r[k])
    return rows

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# what the minimum cluster size (and the other filters) did: the one post-filter result
import anndata
final = anndata.read_h5ad("${h5ad}", backed="r")
after = final.obs["${meta.author_cell_type}"].astype(str).value_counts()
before = {r["cluster"]: int(r["count"]) for r in csv.DictReader(open("${sizes_before}"))}
dropped = [{"cluster": c, "cells": n} for c, n in sorted(before.items()) if c not in after.index]
size = next(csv.DictReader(open("${summary_before}")))
fc = record.get("filter_choices", {})

out = {
    "schema_version": "1.0",
    "kind": "master_dataset_summary",
    "dataset_id": "${meta.dataset_id}",
    "dataset_version_id": "${meta.dataset_version_id}",
    "harvester": record,
    "sc_nsforest_qc": {
        **run,
        "input": {
            "cellxgene_harvester_version": fc.get("harvester_version"),
            "cellxgene_harvester_snapshot": fc.get("snapshot"),
            "harvester_csv": {"file": os.path.basename("${params.harvester_csv}"), "sha256": sha256("input/harvester.csv")},
            "harvester_json": {"file": "${harvester_json}", "sha256": sha256("${harvester_json}")},
            "filtered_h5ad_url": record.get("filtered_h5ad_url"),
            "input_h5ad": json.load(open("${input_info}")),
        },
        "filtering": {
            "min_cluster_size": run["min_cluster_size"],
            "cells_before": int(size["n_obs"]),
            "genes_before": int(size["n_vars"]),
            "clusters_before": int(size["n_clusters"]),
            "cells_after": int(final.n_obs),
            "genes_after": int(final.n_vars),
            "clusters_after": int(len(after)),
            "clusters_dropped": dropped,
        },
        "nsforest_image": "${task.container}",
        "gene_annotation": gene_annotation(),
        "dataset_summary": summary,
        "binary_genes": {
            "positive": table("${positive_csv}"),
            "top10": table("${top10_csv}"),
        },
    },
}
name = "master_dataset_summary_" + "${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}".replace("/", "_").replace(" ", "_") + ".json"
json.dump(out, open(name, "w"), indent=2)
PY
    """
}
