/**
 * Download H5AD Module
 *
 * Downloads an h5ad file from CellxGene (https://) or AWS S3 (s3://).
 * Retries up to 3 times on failure.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with first_author, year (and all other dataset metadata)
 *   - url:  https:// or s3:// URL to the h5ad file
 *
 * Output:
 * -------
 * @emit h5ad: tuple(meta, {first_author}_{year}.h5ad)
 */
process download_h5ad_process {
    tag "download_${meta.first_author}_${meta.year}"

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
