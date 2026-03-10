/**
 * Prep Medians and Binary Scores Module
 *
 * Runs ns.pp.prep_medians() then ns.pp.prep_binary_scores() on the filtered
 * AnnData in a single step, matching DEMO Section 3 exactly.
 * varm entries are extracted to CSV before writing h5ad to avoid h5py
 * TypeError when cluster names contain '/' (e.g. 'VSMC/P').
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta: Map with organ, first_author, year, author_cell_type
 *   - h5ad: Path to adata_filtered.h5ad
 *
 * Output:
 * -------
 * @emit complete: tuple(meta, adata_prep.h5ad, [medians.csv, binary_scores.csv, ...])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_medians.{csv,pkl}
 *                   {organ}_{first_author}_{year}_{cluster_header_safe}_binary_scores.{csv,pkl}
 */
process prep_medians_binary_scores_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_adata_prep.h5ad"),
          path("${meta.organ}_${meta.first_author}_${meta.year}_*.{csv,pkl}", optional: true),
          emit: complete

    script:
    """
    nsforest-cli prep-medians-binary-scores \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
