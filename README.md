>📋  A template README.md for code accompanying a Machine Learning paper

# *_S/M/L-QSD_*: Three quality-associated sequencing datasets to evaluate anomaly detection

This repository is the official implementation of [*_S/M/L-QSD_*: Three quality-associated sequencing datasets to evaluate anomaly detection](Place-holder-url). 

![Workflow to create the datasets](data/qsd_creation_v4.jpg)

## Requirements

To install requirements for the ML Experiments:

```setup
pip install -r experiments/requirements.txt
```

For the Feature generation scirpts: 
```
Docker 
```

## Experimental Results
After running the scripts "unsupervised_experiments.py" and "supervised_experiments.py", you will receive the following performance (AUC ROC mean ± standard deviation) of unsupervised anomaly detection
(top) and supervised classification (bottom) algorithms for the ChIP-Seq assay on S-QSD, M-QSD, and L-QSD.  

![image](https://github.com/user-attachments/assets/8333b844-f5c0-4bd0-a647-ea5754065c9d)


