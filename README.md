
## Nextflow Workflow

**Overview**

This repository includes a modular Nextflow pipeline that executes [NSForest](https://github.com/JCVenterInstitute/NSForest).  This program finds for each cluster specified the necessary and sufficient markers based upon the data as measured and prepared in the single cell RNA seq results in a supplied h5ad file.  Additionally, the top 10 binary genes are also provided.  Each computational step is implemented as an individual Nextflow module in the modules/ directory. The master workflow is defined in main.nf.


## input

The nextflow begins by parsing the input file formated as below:

The contents of this file are a test for testing the workflow, it is the same file, test.csv used in testing the [scsilhouette-nf](https://github.com/nih-nlm/scsilhouette-nf)  Nextflow workflow:

```bash
h5ad_file,label_key,embedding_key,organism,disease,filter_normal,metric,save_scores,save_cluster_summary,save_annotation,tissue,author,publication_date,publication,cell_count
https://datasets.cellxgene.cziscience.com/10b4b990-9f06-441f-be9a-f2dfbd353716.h5ad,cell_type,X_umap,Mus_musculus,normal,True,euclidean,True,True,True,embryo,Sampath_Kumar,2025,Nat_Genet,2567
```

It has a header file that contains:

* h5ad_file location
* label_key - this is the clustering label
* embedding_key - mostly X_umap but ocassionaly not supported so could be X_tsne or other
* organism (e.g. Homo sapiens)
* disease (e.g. normal, neoplasm, etc)
* tissue (e.g. kidney -- actually an organ level specification)
* filter_normal - this is True or False indicating whether or not the file should only be focused on the normal tissues
* metric - not used for this NSForest function - used with scsilhouette score calculation
* save_scores - not used for this NSForest function - used with scsilhouette score calculation
* save_cluster_summary - not used for this NSForest function - used with scsilhouette score calculation
* save_annotation - not used for this NSForest function - used with scsilhouette score calculation
* author (first author of a publication associated with the submitted dataset) - not used
* year (publication year) - not used
* publication (name of publication) - not used
* cell_count (number of cells submitted in the dataset) - not used

## Usage

The execution has been conducted on [Lifebit's free NF-Copilot](cloudos.lifebit.ai)

1. There you select this Nextflow workflow [nsforest-nf](https://github.com/nih-nlm/nsforest-nf)
2. Add your parameter *'dataset_csv'*, upload from your computer or use from the dataset selected with File Explorer
3. Select a low-cost, low memory, low cpu master node as your execution node
4. Select your job queue
5. Select resumable
6. Then Run Analysis.


## Modular Design

Each pipeline process is defined in a separate module in modules/:

```bash
modules/
├── nsforest_process.nf
```


