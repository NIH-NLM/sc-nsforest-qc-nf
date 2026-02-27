#!/usr/bin/env python3
"""
parse_nf_docs.py — Nextflow module documentation generator for Sphinx.

Walks the modules/ directory tree, extracts /** ... */ docblocks and
``process NAME {`` declarations from every .nf file, and writes one
.rst file per module group (nsforest, scsilhouette, publish) into
docs/source/nextflow/.

Run this script before sphinx-build so the generated .rst files are
always in sync with the actual .nf source:

    python3 docs/parse_nf_docs.py
    sphinx-build -b html docs/source docs/build/html

The generated files are intentionally committed to the repository so that
ReadTheDocs and GitHub Pages can build the docs without running this script.
Regenerate them whenever .nf files change.

Usage
-----
From the repository root::

    python3 docs/parse_nf_docs.py [--modules-dir MODULES_DIR]
                                   [--output-dir OUTPUT_DIR]
                                   [--dry-run]

Defaults:
    --modules-dir  modules/
    --output-dir   docs/source/nextflow/
"""

import argparse
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

_DOCBLOCK_RE  = re.compile(r'/\*\*(.*?)\*/', re.DOTALL)
_PROCESS_RE   = re.compile(r'process\s+(\w+)\s*\{')
_PARAM_RE     = re.compile(r'params\.(\w+)', re.MULTILINE)


def _strip_leading_stars(text: str) -> str:
    """Remove the leading ' * ' from each line of a /** ... */ block."""
    lines = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith('* '):
            lines.append(stripped[2:])
        elif stripped == '*':
            lines.append('')
        else:
            lines.append(stripped)
    return '\n'.join(lines).strip()


def parse_nf_file(path: Path) -> dict | None:
    """
    Parse a single .nf file.

    Returns a dict with keys:
        name        str   — process name (e.g. filter_adata_process)
        docstring   str   — cleaned docblock text, or empty string
        params      list  — params.* references found in the file
        source_file str   — relative path used in the .rst for reference
    Returns None if no process declaration is found.
    """
    source = path.read_text(encoding='utf-8')

    # Find process name
    m = _PROCESS_RE.search(source)
    if not m:
        return None
    process_name = m.group(1)

    # Find the last /** ... */ block before the process declaration
    docstring = ''
    for dm in _DOCBLOCK_RE.finditer(source):
        if dm.start() < m.start():
            docstring = _strip_leading_stars(dm.group(1))

    # Collect params references
    params = sorted(set(_PARAM_RE.findall(source)))

    return {
        'name':        process_name,
        'docstring':   docstring,
        'params':      params,
        'source_file': str(path),
    }


# ---------------------------------------------------------------------------
# RST generation helpers
# ---------------------------------------------------------------------------

def _title(text: str, char: str = '-') -> str:
    return f"{text}\n{char * len(text)}\n"


def _code(text: str, language: str = 'nextflow') -> str:
    indent = '\n'.join('   ' + l for l in text.splitlines())
    return f".. code-block:: {language}\n\n{indent}\n"


def _docstring_to_rst(docstring: str) -> str:
    """
    Convert the cleaned docblock text to valid RST.

    Handles:
    - Section headings (lines ending with ':' that are followed by a
      line of '-' or are standalone short labels)
    - Code blocks indented by 4+ spaces → .. code-block:: bash
    - Plain paragraphs
    """
    if not docstring:
        return ''

    lines = docstring.splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Detect heading-like lines: "Input:", "Output:", "Script:", etc.
        if re.match(r'^[A-Z][^.]{0,40}:\s*$', line):
            out.append('')
            out.append(line.rstrip())
            out.append('~' * len(line.rstrip()))
            i += 1
            continue

        # Detect code block: 4+ spaces AND looks like a shell command or path
        # (starts with a command word, path, or bash operator).
        # Avoid treating @param / - bullet continuation lines as code.
        is_code_start = (
            (line.startswith('    ') or line.startswith('\t'))
            and bool(re.match(r'\s+(\$|//|#!|[a-z]+-[a-z]+|/|python|bash|nextflow|git|gh|cp|mkdir)', line))
        )
        if is_code_start:
            code_lines = []
            while i < len(lines) and (lines[i].startswith('    ') or lines[i].startswith('\t') or lines[i].strip() == ''):
                code_lines.append(lines[i])
                i += 1
            while code_lines and code_lines[-1].strip() == '':
                code_lines.pop()
            if code_lines:
                dedented = '\n'.join(l[4:] if l.startswith('    ') else l for l in code_lines)
                out.append('')
                out.append('.. code-block:: bash')
                out.append('')
                for cl in dedented.splitlines():
                    out.append('   ' + cl)
                out.append('')
            continue

        out.append(line)
        i += 1

    return '\n'.join(out).strip()


def process_to_rst(info: dict) -> str:
    """Render a single process as an RST section."""
    parts = []

    # Process name as heading
    label = info['name'].replace('_', ' ').title()
    parts.append(_title(label, '^'))
    parts.append(f".. rubric:: ``{info['name']}``\n")
    parts.append(f"*Source:* ``{info['source_file']}``\n")

    # Docstring body
    rst_doc = _docstring_to_rst(info['docstring'])
    if rst_doc:
        parts.append(rst_doc)
        parts.append('')

    # Params table
    if info['params']:
        parts.append('**Params referenced:**\n')
        for p in info['params']:
            parts.append(f'- ``params.{p}``')
        parts.append('')

    return '\n'.join(parts)


def group_to_rst(group_name: str, processes: list[dict], description: str = '') -> str:
    """Render all processes in one module group as a top-level RST page."""
    title_text = f"{group_name.title()} Modules"
    lines = [
        _title(title_text, '='),
    ]
    if description:
        lines.append(description)
        lines.append('')

    for info in processes:
        lines.append(process_to_rst(info))
        lines.append('')

    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Group descriptions
# ---------------------------------------------------------------------------

GROUP_DESCRIPTIONS = {
    'nsforest': (
        "Nextflow modules for the NSForest marker-gene discovery branch of the\n"
        "``sc-nsforest-qc-nf`` workflow.  These modules run inside the\n"
        "``ghcr.io/nih-nlm/nsforest-nf/nsforest`` container and are\n"
        "orchestrated by ``main.nf``.\n\n"
        "Execution order:\n\n"
        "1. ``filter_adata_process`` — tissue / disease / age ontology filter\n"
        "2. ``dendrogram_process`` — cluster dendrogram\n"
        "3. ``cluster_stats_process`` — drives scatter/gather parallelisation\n"
        "4. ``prep_medians_process`` (parallel per cluster)\n"
        "5. ``merge_medians_process``\n"
        "6. ``prep_binary_scores_process`` (parallel per cluster)\n"
        "7. ``merge_binary_scores_process``\n"
        "8. ``plot_histograms_process``\n"
        "9. ``run_nsforest_process`` (parallel per cluster)\n"
        "10. ``merge_nsforest_results_process``\n"
        "11. ``plots_process``\n"
    ),
    'scsilhouette': (
        "Nextflow modules for the scsilhouette silhouette-score QC branch of the\n"
        "``sc-nsforest-qc-nf`` workflow.  These modules run inside the\n"
        "``ghcr.io/nih-nlm/scsilhouette`` container and are orchestrated\n"
        "by ``main.nf``.\n\n"
        "Execution order (runs in parallel with the NSForest branch):\n\n"
        "1. ``compute_silhouette_process``\n"
        "2. ``viz_summary_process``\n"
        "3. ``viz_distribution_process``\n"
        "4. ``viz_dotplot_process``\n"
        "5. ``compute_summary_stats_process``\n"
    ),
    'publish': (
        "Nextflow module that fires **once**, after every dataset in both the\n"
        "NSForest and scsilhouette branches has completed (or failed).\n\n"
        "It clones ``NIH-NLM/cell-kn``, copies outputs into the\n"
        "``data/prod/{organ}/`` tree, commits, pushes a new branch, and\n"
        "opens a pull request against ``main``.\n\n"
        ".. warning::\n\n"
        "   ``params.github_token`` must be a GitHub personal access token\n"
        "   with ``repo`` write access.  Never hardcode it — pass via\n"
        "   ``-params-file params.json`` or an environment variable.\n"
        "   See the README for the three recommended approaches.\n"
    ),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--modules-dir', default='modules',
                    help='Root directory containing module subdirectories (default: modules/)')
    ap.add_argument('--output-dir', default='docs/source/nextflow',
                    help='Where to write generated .rst files (default: docs/source/nextflow/)')
    ap.add_argument('--dry-run', action='store_true',
                    help='Print what would be written without writing files')
    args = ap.parse_args()

    modules_root = Path(args.modules_dir)
    output_root  = Path(args.output_dir)

    if not modules_root.is_dir():
        print(f"ERROR: modules directory not found: {modules_root}", file=sys.stderr)
        sys.exit(1)

    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)

    # Discover groups (subdirectories of modules/)
    groups = sorted(p for p in modules_root.iterdir() if p.is_dir())
    index_entries = []

    for group_dir in groups:
        group_name = group_dir.name
        nf_files   = sorted(group_dir.glob('*.nf'))

        if not nf_files:
            continue

        processes = []
        for nf_path in nf_files:
            info = parse_nf_file(nf_path)
            if info:
                processes.append(info)

        if not processes:
            continue

        description = GROUP_DESCRIPTIONS.get(group_name, '')
        rst_content = group_to_rst(group_name, processes, description)

        out_path = output_root / f"{group_name}_modules.rst"
        if args.dry_run:
            print(f"[dry-run] Would write {out_path} ({len(rst_content)} chars, "
                  f"{len(processes)} processes)")
        else:
            out_path.write_text(rst_content, encoding='utf-8')
            print(f"Written: {out_path}  ({len(processes)} processes)")

        index_entries.append(f"nextflow/{group_name}_modules")

    # Write the nextflow index page
    nextflow_index = "Nextflow Workflow Modules\n=========================\n\n"
    nextflow_index += (
        "The ``sc-nsforest-qc-nf`` Nextflow workflow is composed of three\n"
        "module groups.  These pages are auto-generated from the\n"
        "``/** ... */`` docblocks in each ``.nf`` file by\n"
        "``docs/parse_nf_docs.py``.\n\n"
        "Re-generate after any ``.nf`` change::\n\n"
        "   python3 docs/parse_nf_docs.py\n\n"
        ".. toctree::\n"
        "   :maxdepth: 2\n\n"
    )
    for entry in index_entries:
        nextflow_index += f"   {entry}\n"

    nf_index_path = output_root / "index.rst"
    if args.dry_run:
        print(f"[dry-run] Would write {nf_index_path}")
    else:
        nf_index_path.write_text(nextflow_index, encoding='utf-8')
        print(f"Written: {nf_index_path}")

    print("Done.")


if __name__ == '__main__':
    main()
