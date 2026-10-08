/**
 * Build the ETL JSON
 *
 * Writes the one JSON file for a dataset that goes on to the ETL. It holds the whole
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
 *   - harvester_json: <dataset_id>.filtered.json from cellxgene-harvester-nf
 *   - h5ad_name:      file name of the final filtered h5ad (adata_filtered_*.h5ad)
 *   - uberon_json, disease_json, hsapdv_json: the resolve files given to the workflow; their
 *     names and sha256 are recorded in sc_nsforest_qc.resolve_files
 *
 * Output:
 * -------
 * @emit json: tuple(meta, sc_nsforest_qc_{prefix}.json)
 */
process build_etl_json_process {
    tag "etl_json_${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'publish_json'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(summary_csvs), path(harvester_json), val(h5ad_name)
    path(uberon_json,  stageAs: 'resolve/uberon.json')
    path(disease_json, stageAs: 'resolve/disease.json')
    path(hsapdv_json,  stageAs: 'resolve/hsapdv.json')

    output:
    tuple val(meta), path("sc_nsforest_qc_*.json"), emit: json

    script:
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
        resolve_file_names:    [uberon: file(params.uberon_json).name, disease: file(params.disease_json).name, hsapdv: file(params.hsapdv_json).name],
        filtered_h5ad:         h5ad_name,
        s3_filtered_h5ad:      params.s3_h5ad_prefix ? "${params.s3_h5ad_prefix.toString().replaceAll('/+\$', '')}/${h5ad_name}" : h5ad_name,
    ])
    """
    cat > run.json <<'RUN_JSON'
${run}
RUN_JSON
python3 - <<'PY'
import csv, glob, hashlib, json
record = json.load(open("${harvester_json}"))
run = json.load(open("run.json"))
names = run.pop("resolve_file_names")
run["resolve_files"] = {
    key: {"file": names[key], "sha256": hashlib.sha256(open("resolve/" + key + ".json", "rb").read()).hexdigest()}
    for key in ("uberon", "disease", "hsapdv")
}
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
out = {
    "schema_version": "1.0",
    "harvester": record,
    "sc_nsforest_qc": {**run, "dataset_summary": summary},
}
name = "sc_nsforest_qc_" + "${meta.organ}_${meta.first_author}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}".replace("/", "_").replace(" ", "_") + ".json"
json.dump(out, open(name, "w"), indent=2)
PY
    """
}
