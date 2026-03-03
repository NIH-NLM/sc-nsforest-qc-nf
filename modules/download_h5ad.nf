process download_h5ad_process {
    tag "${meta.first_author}_${meta.year}"
    label 'download'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'fail' }
    maxRetries 3

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${meta.first_author}_${meta.year}.h5ad"), emit: h5ad

    script:
    """
    echo "Downloading: ${url}"
    wget --quiet --tries=3 --timeout=300 --retry-connrefused \
        -O "${meta.first_author}_${meta.year}.h5ad" \
        "${url}"
    """
}
