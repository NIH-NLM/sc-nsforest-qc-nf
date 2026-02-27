scsilhouette Documentation
==========================

Silhouette score analysis for single-cell clustering quality control.

`scsilhouette` computes silhouette scores to assess clustering quality in
single-cell RNA-seq data and provides integrated visualizations with NSForest
marker discovery results.  It is used as the scsilhouette component of the
`sc-nsforest-qc-nf <https://github.com/NIH-NLM/sc-nsforest-qc-nf>`_
Nextflow workflow.

.. toctree::
   :maxdepth: 2
   :caption: Python API

   modules

.. toctree::
   :maxdepth: 2
   :caption: Nextflow Workflow

   nextflow/index

.. toctree::
   :maxdepth: 1
   :caption: More

   changelog

Installation
------------

.. code-block:: bash

   pip install scsilhouette

Quick Start — standalone
------------------------

.. code-block:: bash

   scsilhouette compute-silhouette \
       --h5ad-path data.h5ad \
       --cluster-header "cell_type" \
       --embedding-key "X_umap" \
       --organ "kidney" \
       --first-author "Lake" \
       --year "2023"

Quick Start — with ontology filtering
--------------------------------------

Generate the three JSON files once using
`cellxgene-harvester <https://github.com/NIH-NLM/cellxgene-harvester>`_:

.. code-block:: bash

   cellxgene-harvester resolve-uberon  kidney  > data/uberon_kidney.json
   cellxgene-harvester resolve-disease normal  > data/disease_normal.json
   cellxgene-harvester resolve-hsapdv  --min-age 15 > data/hsapdv_adult_15.json

Then run with filtering:

.. code-block:: bash

   scsilhouette compute-silhouette \
       --h5ad-path data.h5ad \
       --cluster-header "cell_type" \
       --embedding-key "X_umap" \
       --organ "kidney" \
       --first-author "Lake" \
       --year "2023" \
       --filter-normal \
       --uberon  data/uberon_kidney.json \
       --disease data/disease_normal.json \
       --hsapdv  data/hsapdv_adult_15.json

Quick Start — Nextflow workflow
--------------------------------

.. code-block:: bash

   nextflow run main.nf \
       --datasets_csv data/homo_sapiens_kidney_harvester_final.csv \
       --organ        kidney \
       --uberon_json  data/uberon_kidney.json \
       --disease_json data/disease_normal.json \
       --hsapdv_json  data/hsapdv_adult_15.json \
       --github_token "$(cat ~/.github_token)" \
       -c             configs/macamd64.config

.. warning::

   Pass ``--github_token`` via ``-params-file params.json`` or an environment
   variable — never hardcode it in a config file or on the command line where
   it may appear in shell history.  See the README for full details.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
