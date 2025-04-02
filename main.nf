// main.nf - Nextflow pipeline to run NSForest using the provided Docker image

params.cluster = ''
params.input_h5ad = ''
params.outdir = 'results'

process RunNSForest {
    tag "cluster: ${params.cluster}"
    publishDir "${params.outdir}", mode: 'copy'

    container 'ghcr.io/nih-nlm/nsforest-docker:latest'

    input:
    path input_h5ad from file(params.input_h5ad)

    output:
    path "${params.cluster}_output", emit: output_dir

    script:
    """
    mkdir -p ${params.cluster}_output
    nsforest --input \$input_h5ad --cluster ${params.cluster} --output ${params.cluster}_output
    """
}

workflow {
    RunNSForest()
}

