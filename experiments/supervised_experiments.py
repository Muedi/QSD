import os
import utils
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier

os.environ["CUDA_VISIBLE_DEVICES"] = "-1" # tensorflow only sees CPU (disable use of GPU)
import tensorflow as tf
from joblib import Parallel, delayed


tf.config.set_visible_devices([], "GPU") # disable access to GPU

# load data
file_name = "s_qsd.csv"
data = pd.read_csv(file_name, header = 0, sep = ";")

# remove metadata and convert target column to numeric representation (1: revoked/ 0: released)
preprocessed_data = data.copy()
preprocessed_data = preprocessed_data.loc[preprocessed_data["assay"] == "ChIP-seq"].reset_index(drop = True) # filter assay type
preprocessed_data = preprocessed_data.drop(columns = ["accession", "assay", "organism"], axis = 1) # removes metadata
preprocessed_data["status"] = preprocessed_data["status"].replace("revoked", 1.0) # anomalies
preprocessed_data["status"] = preprocessed_data["status"].replace("released", 0.0) # normal data

# separate features and target column
features = preprocessed_data.copy()
target = preprocessed_data.copy()

features = features.drop(columns = ["status"])
target = preprocessed_data["status"]

# define supervised learning models
models = ["logreg", "rf", "gb", "ann"]

# train and evaluate models in parallel
auc_results = Parallel(n_jobs=30)(
    delayed(utils.worker_supervised_experiments)(i, model_name, features, target) 
    for model_name in models
    for i in range(10)) # number of runs

# save results as a DataFrame
df_auc_results = pd.DataFrame(auc_results)

# present results a a latex table
df_summary_supervised_results = (df_auc_results.groupby("Model")[["AUC_ROC_Train", "AUC_ROC_Test",
                                                                   "AUC_PR_Train", "AUC_PR_Test"]]).agg(["mean", "std"])

df_formatted_supervised_results = pd.DataFrame({
    "Model": df_summary_supervised_results.index,
    "AUC_ROC_Train": utils.format_mean_std(df_summary_supervised_results["AUC_ROC_Train"]),
    "AUC_ROC_Test": utils.format_mean_std(df_summary_supervised_results["AUC_ROC_Test"]),
    "AUC_PR_Train": utils.format_mean_std(df_summary_supervised_results["AUC_PR_Train"]),
    "AUC_PR_Test": utils.format_mean_std(df_summary_supervised_results["AUC_PR_Test"]),
})

latex_table_supervised_results = df_formatted_supervised_results.to_latex(index = False, escape = False)
print(latex_table_supervised_results)







