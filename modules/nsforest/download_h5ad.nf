process download_h5ad_process {
    tag "${meta.first_author}_${meta.year}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'fail' }
    maxRetries 3

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${meta.first_author}_${meta.year}.h5ad"), emit: h5ad

    script:
    """
    echo "Downloading: ${url}"
    if [[ "${url}" == s3://* ]]; then
        aws s3 cp "${url}" "${meta.first_author}_${meta.year}.h5ad"
    else
        curl -fL "${url}" > "${meta.first_author}_${meta.year}.h5ad"
    fi
    """
}
