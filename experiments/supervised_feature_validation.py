import os
import utils
import pandas as pd

os.environ["CUDA_VISIBLE_DEVICES"] = "-1" # tensorflow only sees CPU (disable use of GPU)
import tensorflow as tf
from joblib import Parallel, delayed

tf.config.set_visible_devices([], "GPU") # disable access to GPU

# load data
file_name = "S_QSD.csv"
data = pd.read_csv(file_name, header = 0, sep = ",")

df_train_test = pd.read_csv("by_assay_splits/splits_ChIP-seq.csv")

df_train_test = df_train_test.rename(columns = {"Accession": "accession"})
df_train_test = df_train_test.drop("Dataset", axis = 1)


# add columns based on accesssion
df_merge = pd.merge(df_train_test,data, how = "inner", on = "accession" )

# remove metadata and convert target column to numeric representation (1: revoked/ 0: released)
preprocessed_data = df_merge.copy()
preprocessed_data = preprocessed_data.loc[preprocessed_data["assay"] == "ChIP-seq"].reset_index(drop= True) # filter assay type
preprocessed_data = preprocessed_data.drop(columns = ["accession", "assay", "organism"], axis = 1) # removes metadata
preprocessed_data["status"] = preprocessed_data["status"].replace("revoked", 1.0) # anomalies
preprocessed_data["status"] = preprocessed_data["status"].replace("released", 0.0) # normal data

# separate features and target column
features = preprocessed_data.copy()
target = preprocessed_data.copy()

# train and evaluate models in parallel
models = ["logreg", "logreg_balanced", "rf", "rf_balanced", "gb", "gb_balanced", "ann", "ann_focal", "ann_balanced_keras"]
splits = ["split_1", "split_2", "split_3", "split_4", "split_5", "split_6", "split_7", "split_8", "split_9", "split_10"]

# train and evaluate models in parallel
auc_results = Parallel(n_jobs=1)(
    delayed(utils.worker_supervised_experiments)(i, splits, model_name, features, target) 
    for model_name in models
    for i in range(10)) # number of runs

# save results as a DataFrame
df_auc_results = pd.DataFrame(auc_results)

# present results as a latex table 
custom_order = ["logreg", "logreg_balanced", "rf", "rf_balanced", "gb", "gb_balanced", "ann", "ann_focal", "ann_balanced_keras"]

df_summary_supervised_results = (df_auc_results.groupby("Model")[["AUC_ROC_Train", "AUC_ROC_Test",
                                                                   "AUC_PR_Train", "AUC_PR_Test"]]).agg(["mean", "std"]).reindex(custom_order)

df_formatted_supervised_results = pd.DataFrame({
    "Model": df_summary_supervised_results.index,
    "AUC_ROC_Train": utils.format_mean_std(df_summary_supervised_results["AUC_ROC_Train"]),
    "AUC_ROC_Test": utils.format_mean_std(df_summary_supervised_results["AUC_ROC_Test"]),
    "AUC_PR_Train": utils.format_mean_std(df_summary_supervised_results["AUC_PR_Train"]),
    "AUC_PR_Test": utils.format_mean_std(df_summary_supervised_results["AUC_PR_Test"]),
})


latex_table_supervised_results = df_formatted_supervised_results.to_latex(index = False, escape = False)
print(latex_table_supervised_results)

