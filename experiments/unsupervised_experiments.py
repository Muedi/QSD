import torch
import utils
import pandas as pd
import tensorflow as tf
from pyod.models.knn import KNN
from pyod.models.lof import LOF
from pyod.models.iforest import IForest
from pyod.models.auto_encoder import AutoEncoder
from pyod.models.vae import VAE
from pyod.models.deep_svdd import DeepSVDD
from joblib import Parallel, delayed


tf.config.set_visible_devices([], "GPU") # disable access to GPU

# load data
file_name = "s_qsd.csv"
data = pd.read_csv(file_name, header = 0, sep = ";")

# remove metadata and convert target column to numeric representation (1: revoked/ 0: released)
preprocessed_data = data.copy()
preprocessed_data = preprocessed_data.loc[preprocessed_data["assay"] == "ChIP-seq"].reset_index(drop= True) # filter assay type
preprocessed_data = preprocessed_data.drop(columns = ["accession", "assay", "organism"], axis = 1) # removes metadata
preprocessed_data["status"] = preprocessed_data["status"].replace("revoked", 1.0) # anomalies
preprocessed_data["status"] = preprocessed_data["status"].replace("released", 0.0) # normal data

# separate features and target column
features = preprocessed_data.copy()
target = preprocessed_data.copy()

features = features.drop(columns = ["status"], axis = 1)
target = preprocessed_data["status"]

# define unsupervised learning models
models = [KNN(),
          LOF(), 
          IForest(),
          AutoEncoder(device ="cpu", random_state=42),
          VAE(device ="cpu", random_state=42),
          DeepSVDD(n_features = len(features.columns), random_state=42)]

# train and evaluate models in parallel
auc_results = Parallel(n_jobs=30)(
    delayed(utils.worker_unsupervised_experiments)(i, model, features, target)
    for model in models
    for i in range(10) # number of runs
)

# save results as a DataFrame
df_auc_results = pd.DataFrame(auc_results) 

# present results as a latex table 
df_summary_unsupervised_results = (df_auc_results.groupby("Model")[["AUC_ROC", "AUC_PR"]]).agg(["mean", "std"])

df_formatted_unsupervised_results = pd.DataFrame({
    "Model": df_summary_unsupervised_results.index,
    "AUC_ROC": utils.format_mean_std(df_summary_unsupervised_results["AUC_ROC"]),
    "AUC_PR": utils.format_mean_std(df_summary_unsupervised_results["AUC_PR"]),
})

latex_table_unsupervised_results = df_formatted_unsupervised_results.to_latex(index = False, escape = False)
print(latex_table_unsupervised_results)





