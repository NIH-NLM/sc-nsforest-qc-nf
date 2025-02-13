#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage: nextflow rnu nsforest-parallel-by-cluster.nf --csvPath '../data/cellxgene-sample/cellxgene-sample.csv'

    Options:
    --help
        Show help message
    """.stripIndent()
}

params.csvPath = ""
params.help = ""

// Show help message, if requested
if (params.help) {
  helpMessage()
  exit 0
}

/* Create tuples of input data files and corresponding cluster headers
for preprocessing */
preprocess_input_tuple_ch = channel
    .fromPath(params.csvPath)
    .splitText()
    .splitCsv()

// Define the expected NS-Forest output CSV file header
results_header_ch = channel.value(
    "clusterName,clusterSize,f_score,PPV,recall,TN,FP,FN,TP,marker_count,NSForest_markers,binary_genes,onTarget"
)

process preprocess_cxg_file {

    publishDir "results", mode: "copy"

    input:
    tuple path(cxg_file), val(cluster_header)
    val results_header

    output:
    tuple env(key), path("*_pp.*"), val(cluster_header), path("${cxg_file.baseName}_${cluster_header}_results.csv")
    path "clusters.txt"
    path "${cxg_file.baseName}_${cluster_header}_results.csv"
    path "${cxg_file.baseName}_pp.${cxg_file.extension}"

    script:
    """
    # Preprocess the specified file
    nsforest.py --preprocess-adata-file -c "${cluster_header}" "${cxg_file}"

    # Write header to results file
    echo "${results_header}" > "${cxg_file.baseName}_${cluster_header}_results.csv"

    # Create key for join
    key=\$(tail -1 clusters.txt | cut -d "," -f 1)
    """

}

process run_nsforest_by_cluster {

    publishDir "results", mode: "copy"

    input:
    tuple val(key), path(pp_cxg_file), val(cluster_header), path(collected_results_file), val(cluster_name)
    val results_header

    output:
    path collected_results_file

    script:
    """
    # Run NSforest by cluster using the preprocessed file
    nsforest.py --run-nsforest-without-preprocessing -c "${cluster_header}" -l "${cluster_name}" "${pp_cxg_file}"

    # Collect the cluster results ignoring the results header
    cat "${cluster_header}_results.csv" | grep -v "${results_header}" >> "${collected_results_file}"
    """
}

workflow {

    /* Preprocess the input file, noting the files that will be used
    to collect the results by cluster */
    (
        preprocess_output_tuple_ch, cluster_names_files_ch, collected_results_files_ch, preprocessed_files_ch
    ) = preprocess_cxg_file(
        preprocess_input_tuple_ch, results_header_ch
    )

    // Create tuples of preprocessed data files, cluster header, and output data files
    cluster_names_ch = cluster_names_files_ch
        .splitText()
        .splitCsv()
    nsforest_input_tuple_ch = preprocess_output_tuple_ch
        .combine(cluster_names_ch, by: 0)

    /* Run NSForest on each input file and cluster header, for each
    cluster */
    run_nsforest_by_cluster(nsforest_input_tuple_ch, results_header_ch)

}
