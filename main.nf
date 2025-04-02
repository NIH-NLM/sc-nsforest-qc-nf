// main.nf - Final Safe Version for NSForest with spaces-safe workaround

nextflow.enable.dsl=2

params.cluster = ''
params.input_h5ad = ''
params.outdir = 'results'

workflow {
    input_h5ad_file = file(params.input_h5ad)
    nsforest_out = RunNSForest(input_h5ad_file)
    nsforest_out.view()
}

process RunNSForest {

    tag "cluster: ${params.cluster}"
    publishDir "${params.outdir}", mode: 'copy'
    container 'ralatsdio/nsforest:latest'

    input:
    path input_h5ad

    output:
    path "${params.cluster}_output"

    script:
    """
    # Copy input to path without spaces
    safe_input="input.h5ad"
    cp "${input_h5ad}" "\$safe_input"

    mkdir -p "${params.cluster}_output"
    nsforest.py --preprocess-adata-file -c "${params.cluster}" "\$safe_input"
    """
}

