// File: modules/log_nsforest_version.nf

process log_nsforest_version_process {

    tag "log-nsforest-version"

    output:
    path "nsforest_version.txt", emit: nsforest_version

    container:
    'ghcr.io/nih-nlm/nsforest-nf/nsforest:latest'

    publishDir "${params.outdir}/metadata", pattern: "nsforest_version.txt", mode: 'copy'

    script:
    """
    cd /opt/nsforest
    git rev-parse HEAD > nsforest_version.txt || echo 'unknown' > nsforest_version.txt
    """
}
