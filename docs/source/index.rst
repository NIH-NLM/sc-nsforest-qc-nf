sc-nsforest-qc-nf Documentation
================================

Nextflow pipeline for NSForest marker gene discovery and silhouette score
quality control of single-cell RNA-seq data.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   README

Workflow Overview
-----------------

``sc-nsforest-qc-nf`` runs two analysis packages in sequence:

- `NSForest <https://github.com/JCVenterInstitute/NSForest>`_ — marker gene
  discovery by the J. Craig Venter Institute. NSForest is called via
  ``nsforest-cli``, a workflow-internal command-line wrapper bundled inside
  the ``ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest`` container.

- `scsilhouette <https://nih-nlm.github.io/scsilhouette/>`_ — silhouette score
  quality control developed by NIH-NLM. See the
  `scsilhouette documentation <https://nih-nlm.github.io/scsilhouette/>`_ for
  full API reference.

Quick Start
-----------

.. code-block:: bash

   nextflow run main.nf \
       --datasets_csv homo_sapiens_kidney_harvester_final.csv \
       --organ kidney \
       --uberon_json uberon_kidney.json \
       --outdir results/kidney

Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``--datasets_csv``
     - required
     - Path to ``homo_sapiens_{organ}_harvester_final.csv``
   * - ``--organ``
     - required
     - Organ label (e.g. ``kidney``, ``heart``) — used in output directory naming
   * - ``--uberon_json``
     - required
     - Path to ``uberon_{organ}.json`` from cellxgene-harvester resolve-uberon
   * - ``--min_age``
     - ``15``
     - Minimum donor age in years for adult cell filtering
   * - ``--min_cluster_size``
     - ``5``
     - Minimum cells per cluster — smaller clusters dropped and logged
   * - ``--outdir``
     - ``./results``
     - Output directory
   * - ``--publish_mode``
     - ``copy``
     - Nextflow ``publishDir`` mode

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
