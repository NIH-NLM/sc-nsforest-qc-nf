sc-nsforest-qc-nf Documentation
================================

``sc-nsforest-qc-nf`` is a Nextflow pipeline for NSForest marker gene
discovery and silhouette score quality control of single-cell RNA-seq data.

It orchestrates parallel execution of
`NSForest <https://github.com/JCVenterInstitute/NSForest>`_
marker discovery and
`scsilhouette <https://github.com/NIH-NLM/scsilhouette>`_
clustering quality control across multiple datasets and organs, with
ontology-based cell filtering driven by
`cellxgene-harvester <https://github.com/NIH-NLM/cellxgene-harvester>`_
outputs.

It is part of the `NIH NLM Cell Knowledge Network <https://github.com/NIH-NLM/cell-kn>`_.

.. toctree::
   :maxdepth: 2
   :caption: Overview

   README

.. toctree::
   :maxdepth: 2
   :caption: Nextflow Workflow

   nextflow/index

.. toctree::
   :maxdepth: 2
   :caption: Python API (nsforest-cli)

   modules

Related Projects
----------------

- `NSForest <https://github.com/JCVenterInstitute/NSForest>`_ — marker gene discovery algorithm (J. Craig Venter Institute). `NSForest documentation <https://nsforest.readthedocs.io>`_
- `scsilhouette <https://github.com/NIH-NLM/scsilhouette>`_ — silhouette score QC package (NIH NLM). `scsilhouette documentation <https://nih-nlm.github.io/scsilhouette>`_
- `cellxgene-harvester <https://github.com/NIH-NLM/cellxgene-harvester>`_ — single-cell data aggregation from CellxGene
- `cell-kn <https://github.com/NIH-NLM/cell-kn>`_ — NIH NLM Cell Knowledge Network

Quick Start — Nextflow workflow
--------------------------------

.. code-block:: bash

   nextflow run main.nf \
       --datasets_csv data/homo_sapiens_kidney_harvester_final.csv \
       --organ        kidney \
       --uberon_json  data/uberon_kidney.json \
       --disease_json data/disease_normal.json \
       --hsapdv_json  data/hsapdv_adult_15.json \
       --outdir       results/kidney \
       -c             configs/macamd64.config

.. warning::

   Pass ``--github_token`` via ``-params-file params.json`` or an environment
   variable — never hardcode it in a config file or on the command line where
   it may appear in shell history.  See the README for full details.

Repository Structure
--------------------

.. code-block:: text

   sc-nsforest-qc-nf/
   ├── configs/                        # Platform-specific Nextflow configs
   │   ├── aws.config
   │   ├── macamd64.config
   │   └── nexflow_biowulf.config
   ├── container/nsforest/             # Docker image for nsforest-cli
   │   ├── Dockerfile
   │   └── context/src/nsforest_cli/   # nsforest-cli Python package
   ├── docs/                           # Sphinx documentation
   │   ├── parse_nf_docs.py            # Auto-generates RST from .nf docblocks
   │   └── source/
   ├── modules/
   │   ├── nsforest/                   # NSForest Nextflow process modules
   │   ├── publish/                    # cell-kn publish module
   │   └── scsilhouette/               # scsilhouette Nextflow process modules
   ├── main.nf                         # Pipeline entry point
   └── nextflow.config                 # Default parameters and container config

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
