# *_S/M/L-QSD_*: Three quality-associated sequencing datasets to evaluate anomaly detection

This repository is the official implementation of [*_S/M/L-QSD_*: Three quality-associated sequencing datasets to evaluate anomaly detection](Place-holder-url). 

![Workflow to create the datasets](data/qsd_creation_v13.jpg)
We generate feature sets that contain quality information to predict low-quality in fucntional genomics files, we provide a dataset of 37,491 RNA-seq, DNAse-seq, ChIP-seq and eCLIP files. 
We test Anomaliy-detection algorithms and compare them to supervised ML performance of identifying low quality features, based on ENCODE's *released* and *revoked* status. 
## Requirements

#### To install requirements for the ML Experiments:

```setup
pip install -r experiments/requirements.txt
```

#### To install requirements for the feature generation:
```setup
conda env create -f bioinf.yml
```
By using the provided yml file, it is possible to create a conda environment with all dependencies and Bioconductor packages needed. The yml file also provides an overview of all dependencies with their versions.

It can be advisable to follow the more detailed step-by-step guides provided in the [seqQscorer README](https://github.com/salbrec/seqQscorer).

The easiest and fastest way to get ready for using the feature generation pipeline is pulling the seqQscorer docker image, running it and using the feature generation script within the docker.
For more details about how to use the docker image under Linux and Windows systems, we refer to the [seqQscorer README](https://github.com/salbrec/seqQscorer).

## Experimental Results
After running the scripts "unsupervised_feature_validation.py" and "supervised_feature_validation.py", you will receive the following performance (AUC ROC mean ± standard deviation) of unsupervised anomaly detection
(top) and supervised classification (bottom) algorithms for the ChIP-Seq assay on QC-34, and BL-n.  

<img width="601" height="236" alt="image" src="https://github.com/user-attachments/assets/d4fe4b45-6f8b-4bfa-9a13-a4e8379b62af" />
<img width="639" height="228" alt="image" src="https://github.com/user-attachments/assets/27e811ec-6e75-4f46-8dd9-2bd5068ee08d" />


## Feature Generation
The script `feature_generation_pipeline_SML_QSD.py` has only one command line argument for the ENCODE accession of a raw sequencing file in FASTQ format. The pipeline autmoatically downloads this file into `./data/fastq/` and runs different tools and scripts to derive the different features. The result files from this will be saved under different directories in `./data/features/`. This is an example for creating all features for the QC-34, and BL-n datasets for the sample ENCFF001NAO from ENCODE: 
```
python feature_generation_pipeline_SML_QSD.py ENCFF001NAO
```
The feature generation might be applied to multiple samples. To create the dataset containing the features for multiple NGS samples, the `gather_features_for_S_M_L_QSD_datasets.py` script can be used. It gathers all information from the different feature sets to assemble the QC-34, and BL-n datasets and create a comma-separated file for each of them.

## BL-n Generation

The required 'blocklist_read_counts.csv' can be downloaded here: https://seafile.rlp.net/seafhttp/f/7a7db59b51674ea8b7eb/?op=view 









