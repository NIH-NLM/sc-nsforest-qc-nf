import os
import sys

project = 'sc-nsforest-qc-nf'
copyright = '2026, NIH-NLM'
author = 'NIH-NLM'
release = '1.0.0'

extensions = [
    'myst_parser',
    'sphinx.ext.autosectionlabel',
]

# MyST parser settings — enable anchors for cross-referencing headings
myst_heading_anchors = 3

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Make the page title match the repo name rather than the README h1
html_title = 'sc-nsforest-qc-nf'
