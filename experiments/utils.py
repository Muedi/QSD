import os
import time
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, auc, precision_recall_curve
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras import layers
from tensorflow.keras.utils import set_random_seed
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.utils import class_weight


def create_ann(input_dimension):
    """
    This function describes the archictecture and the compiling process of our artifical neural network. 
    
    Parameters
    -----------------
    input_dimension: number of features
    """

    set_random_seed(812) # enhances reproducibility of results

    model = tf.keras.models.Sequential()
    model.add(layers.Dense(50, activation = "relu", input_dim = input_dimension, name = "InputLayer_HiddenLayer0")) 
    model.add(layers.Dense(20, activation = "relu", name = "HiddenLayer1")) 
    model.add(layers.Dense(10, activation = "relu", name = "HiddenLayer2")) 
    model.add(layers.Dense(5, activation = "relu", name = "HiddenLayer3")) 
    model.add(layers.Dense(1, activation = "sigmoid", name = "OutputLayer")) 

    model.compile(loss = "binary_crossentropy",
                optimizer = "adam",
                metrics= [tf.keras.metrics.AUC(curve = "ROC"), tf.keras.metrics.AUC(curve = "PR")])
    return model


def impute_train_test_features_supervised(df_X_train_scaled, df_X_test_scaled, file_name):
   """
   This function imputes or replaces the missing values of the training and testing data of the supervised feature validation. 
   The imputation method depends on the type of the datasets.

   *QC-34*: 
   All missing values are imputed using the median of the training column, except the columns that start with "LOC_". Those "LOC"-features are filled with a constant value (default 0.0).

   *BL-n*: 
   All missing values are imputed using the median of the training column.

  The impuation methods are first fitted on the features of the training data and then applied to the features of the testing data.
   
   
  Parameters
  -----------------
  df_X_train_scaled:  scaled training data
  df_X_test_scaled: scaled testing data 
  file_name: name of the dataset

  """
   is_qc_file = "QC" in file_name

   if is_qc_file:
    orig_train_cols = df_X_train_scaled.columns
    orig_test_cols = df_X_test_scaled.columns
    qc_loc_features = [loc_column for loc_column in df_X_train_scaled if loc_column.startswith("LOC_")]
    qc_other_features = [other_column for other_column in df_X_train_scaled if not other_column.startswith("LOC_")]

    imp_loc  = SimpleImputer(strategy="constant", fill_value=0)
    imp_other = SimpleImputer(strategy="median")

    X_train_loc   = imp_loc.fit_transform(df_X_train_scaled[qc_loc_features])
    X_train_other  = imp_other.fit_transform(df_X_train_scaled[qc_other_features])

    X_test_loc   = imp_loc.transform(df_X_test_scaled[qc_loc_features])
    X_test_other  = imp_other.transform(df_X_test_scaled[qc_other_features])

    df_X_train_scaled_imputed = df_X_train_scaled.copy()
    df_X_test_scaled_imputed = df_X_test_scaled.copy()

    df_X_train_scaled_imputed[qc_loc_features]  = X_train_loc
    df_X_train_scaled_imputed[qc_other_features] = X_train_other

    df_X_test_scaled_imputed[qc_loc_features]  = X_test_loc
    df_X_test_scaled_imputed[qc_other_features] = X_test_other

    assert df_X_train_scaled_imputed.columns.equals(orig_train_cols)
    assert df_X_test_scaled_imputed.columns.equals(orig_test_cols)

   else:
    imp = SimpleImputer(missing_values=np.nan, strategy='median')
    X_train_imputed =  imp.fit_transform(df_X_train_scaled)
    X_test_imputed = imp.transform(df_X_test_scaled)

    df_X_train_scaled_imputed = pd.DataFrame(X_train_imputed, columns = df_X_train_scaled.columns, index = df_X_train_scaled.index)
    df_X_test_scaled_imputed = pd.DataFrame(X_test_imputed, columns = df_X_test_scaled.columns, index = df_X_test_scaled.index)

   return df_X_train_scaled_imputed, df_X_test_scaled_imputed

def impute_features_unsupervised(df_X_scaled, file_name):
   """
   This function imputes or replaces the missing values of the data of the unsupervised feature validation. 
   The imputation method depends on the type of the datasets.

   *QC-34*: 
   All missing values are imputed using the median of the training column, except the columns that start with "LOC_". Those "LOC"-features are filled with a constant value (default 0.0).

   *BL-n*: 
   All missing values are imputed using the median of the training column.
   
   
  Parameters
  -----------------
  df_X_scaled:  scaled input features
  file_name: name of the dataset

  """
   is_qc_file = "QC" in file_name

   if is_qc_file:
    orig_cols = df_X_scaled.columns
    qc_loc_features = [loc_column for loc_column in df_X_scaled if loc_column.startswith("LOC_")]
    qc_other_features = [other_column for other_column in df_X_scaled if not other_column.startswith("LOC_")]

    imp_loc  = SimpleImputer(strategy="constant", fill_value=0)
    imp_other = SimpleImputer(strategy="median")

    X_loc   = imp_loc.fit_transform(df_X_scaled[qc_loc_features])
    X_other  = imp_other.fit_transform(df_X_scaled[qc_other_features])

    df_X_scaled_imputed = df_X_scaled.copy()
    
    df_X_scaled_imputed[qc_loc_features]  = X_loc
    df_X_scaled_imputed[qc_other_features] = X_other

    assert df_X_scaled_imputed.columns.equals(orig_cols)

   else:
    imp = SimpleImputer(missing_values=np.nan, strategy='median')
    X_imputed =  imp.fit_transform(df_X_scaled)
    df_X_scaled_imputed = pd.DataFrame(X_imputed, columns = df_X_scaled.columns, index = df_X_scaled.index)

   return df_X_scaled_imputed


def scale_and_impute_features_supervised_experiments(X_train, X_test, y_train, y_test, file_name):
  """
  This function first checks if there are missing values in any of the columns. 
  If there are none, the features of the training and testing data are transformed or scaled directly with min-max-scaling. 
  The scaler is fitted on the features of the training data and then applied to the features of the testing data.
  If there are missing values, then each column of the training and testing data is first transformed with min-max-scaling before 
  the missing values are replaced or imputed using the function "impute_train_test_features". 
  Again the scaling and impuation methods are first fitted on the features of the training data and then applied to the features of the testing data.
  
  Parameters
  -----------------
  X_train: input features of the training data
  X_test: input features of the testing data
  y_train: target variable of the training data
  y_test: target variable of the testing data
  """
  
  scaler = MinMaxScaler()
  missing_values = X_train[X_train.isna().any(axis = 1)] # filters dataframe to select any row with missing values

# there are no missing values
  if missing_values.shape[0] == 0: 
    
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    df_X_train_scaled_imputed = pd.DataFrame(X_train_scaled, columns = X_train.columns, index = X_train.index)
    df_X_test_scaled_imputed= pd.DataFrame(X_test_scaled, columns = X_test.columns, index = X_test.index)

    df_X_train_scaled_imputed = df_X_train_scaled_imputed.sort_index()
    df_X_test_scaled_imputed = df_X_test_scaled_imputed.sort_index()

  # there are missing values: scale X_train and X_test based on X_train
  else:
    X_train_scaled = scaler.fit_transform(X_train) # MinMaxScaler ignores missing values
    X_test_scaled = scaler.transform(X_test) 
    df_X_train_scaled = pd.DataFrame(X_train_scaled, columns = X_train.columns, index = X_train.index)
    df_X_test_scaled = pd.DataFrame(X_test_scaled, columns = X_test.columns, index = X_test.index)

    # impute missing values with median for each column (imputation of X_train and X_test based on X_train), LOC features for QC-34 need imputation with 0
    df_X_train_scaled_imputed, df_X_test_scaled_imputed = impute_train_test_features_supervised(df_X_train_scaled, df_X_test_scaled, file_name)

    df_X_train_scaled_imputed = df_X_train_scaled_imputed.sort_index()
    df_X_test_scaled_imputed = df_X_test_scaled_imputed.sort_index()

    # add the target column back to the scaled features to ensure that the target and input data are in the same order
    df_X_train_scaled_imputed["status"] = y_train
    df_X_test_scaled_imputed["status"] = y_test

    # test if the order of the entries in the target column is the same after scaling and imputation
    if df_X_train_scaled_imputed["status"].equals(y_train.sort_index()) & df_X_test_scaled_imputed ["status"].equals(y_test.sort_index()):
       df_X_train_scaled_imputed = df_X_train_scaled_imputed.drop(columns = ["status"], axis = 1)
       df_X_test_scaled_imputed = df_X_test_scaled_imputed.drop(columns = ["status"], axis = 1)

    else:
      print("There is an error, check the indices")

  return df_X_train_scaled_imputed, df_X_test_scaled_imputed, y_train.sort_index(), y_test.sort_index()

def scale_and_impute_features_unsupervised_experiments(X, y, file_name):
  """
  This function first checks if there are missing values in any of the columns. 
  If there are none, the features are transformed or scaled directly with min-max-scaling. 
  If there are missing values, then each column is first transformed with min-max-scaling before 
  the missing values are replaced or imputed using the function "impute_train_test_features". 
  
  Parameters
  -----------------
  X: input features
  y: target variable
  """
  scaler = MinMaxScaler()
  df_missing_values = X[X.isna().any(axis = 1)] # filters dataframe to select any row with missing values

 # there are no missing values
  if df_missing_values.shape[0] == 0: 
    X_scaled = scaler.fit_transform(X)
    df_X_scaled_imputed  = pd.DataFrame(X_scaled, columns= X.columns, index=X.index)
    df_X_scaled_imputed = df_X_scaled_imputed.sort_index()
  
  # there are missing values
  else:
    X_scaled = scaler.fit_transform(X) # MinMaxScaler ignores missing values
    df_X_scaled = pd.DataFrame(X_scaled, columns = X.columns, index = X.index)

    # impute missing values with median for each column
    df_X_scaled_imputed = impute_features_unsupervised(df_X_scaled, file_name)

    df_X_scaled_imputed  = df_X_scaled_imputed.sort_index()

    # add the target column back to the scaled features to ensure that the target and input data are in the same order
    df_X_scaled_imputed["status"] = y

    # test if the order of the entries in the target column is the same after scaling and imputation
    if df_X_scaled_imputed["status"].equals(y.sort_index()):
      df_X_scaled_imputed = df_X_scaled_imputed.drop(columns = ["status"], axis = 1)

    else:
      print("There is an error, check the indices")

  return df_X_scaled_imputed , y.sort_index()


def worker_supervised_experiments(i, splits, model_name, features, target, file_name):
  split = splits[i]

  irrelevant_splits = splits.copy() 
  irrelevant_splits.remove(split) # all splits except current split

  # do train test split based on current split
  columns_to_drop = ["status"] + irrelevant_splits
  X = features.drop(columns = columns_to_drop, axis = 1)
  y = target[["status", split]]

  X_train = X.loc[X[split]== 0].reset_index(drop = True)
  X_test = X.loc[X[split]== 1].reset_index(drop = True)

  y_train = y.loc[y[split] == 0].reset_index(drop  = True)
  y_test = y.loc[y[split] == 1].reset_index(drop = True)

  # drop split
  X_train = X_train.drop(columns = [split], axis = 1)
  X_test = X_test.drop(columns = [split], axis = 1)

  y_train = y_train["status"]
  y_test = y_test["status"]

  # re-complete irrelevant splits list 
  irrelevant_splits = splits.copy()

  X_train, X_test, y_train, y_test = scale_and_impute_features_supervised_experiments(X_train, X_test, y_train, y_test, file_name)

  print(f"Split : {split},\n Y-train:{y_train.unique()},\n Y-test{y_test.unique()}\n len x train {len(X_train)}\n len x test {len(X_test)}\n len y train {len(y_train)}\n len y test {len(y_test)}")

  if model_name == "logreg":
     model = LogisticRegression(random_state = 42,max_iter = 1000) 
  elif model_name == "rf":
      model = RandomForestClassifier(random_state = 42)
  elif model_name =="gb":
      model = GradientBoostingClassifier(random_state = 42)
        
  elif model_name == "ann":
      model = create_ann(X_train.shape[1])

  else:
      raise ValueError("Unknown Model")
    
        
  if model_name == "ann":

      start_training = time.time()
      model.fit(X_train, y_train, epochs = 50)
      end_training = time.time()

      y_pred_train = model.predict(X_train).flatten()

      start_inference = time.time()
      y_pred_test= model.predict(X_test).flatten()
      end_inference = time.time()

      training_duration = end_training - start_training # seconds
      inference_duration = end_inference - start_inference # seconds
    
  else:

      start_training = time.time()
      model.fit(X_train, y_train)
      end_training = time.time()

      y_pred_train = model.predict_proba(X_train)[:,1] # outputs the probability that a training sample is classified as "revoked"

      start_inference = time.time()
      y_pred_test = model.predict_proba(X_test)[:,1] # outputs the probability that a testing sample is classified as "revoked"
      end_inference = time.time()

      training_duration = end_training - start_training # seconds
      inference_duration = end_inference - start_inference # seconds

  precision_train, recall_train, _ = precision_recall_curve(y_train, y_pred_train)
  precision_test, recall_test, _ = precision_recall_curve(y_test, y_pred_test)

  return {
        "Run": i,
        "Split": split,
        "Model": model_name,
        "Training_Duration(s)": training_duration, 
        "Inference Time(s)": inference_duration,
        "AUC_ROC_Train": roc_auc_score(y_train, y_pred_train),
        "AUC_ROC_Test": roc_auc_score(y_test, y_pred_test),
        "AUC_PR_Train": auc(recall_train, precision_train),
        "AUC_PR_Test": auc(recall_test, precision_test)
      
  }


def worker_unsupervised_experiments(i, model,  features, target, file_name):
    """
    This function defines the task to be executed via parallel computation. The task covers the 
    scaling and imputation of the features for training the unsupervised learning model, as well as the 
    training and evaluation process of the unsupervised learning model.
    Lastly, the function returns the number of the training and evaluation run, 
    the name of the model and the computed AUC ROC and AUC PR values.

    Parameters
    -----------------
    i: training and evaluation run
    model: unsupervised learning model
    features: input features
    target: target variable
    """
    X = features
    y = target
    
    X, y= scale_and_impute_features_unsupervised_experiments(X,y, file_name)

    start_training = time.time()
    model.fit(X)
    end_training = time.time()

    start_inference = time.time()
    #y_proba = model.predict_proba(X)[:,1] # outputs the probability that a sample is classified as "revoked"
    y_proba = model.decision_scores_ # outputs the score that a sample is classified as "revoked"
    end_inference = time.time()

    training_duration = end_training - start_training # seconds
    inference_duration = end_inference - start_inference # seconds

    precision, recall, _ = precision_recall_curve(y, y_proba)

    return {
        "Run": i, 
        "Model": model.__class__.__name__,
        "Training_Duration(s)": training_duration, 
        "Inference Time(s)": inference_duration,
        "AUC_ROC": roc_auc_score(y, y_proba),
        "AUC_PR": auc(recall, precision)
    }




def format_mean_std(column):
   """
   This function returns the entries in the following format: mean ± standard deviation with 3 decimal places
   """
   return column.apply(lambda row: f"{row['mean']:.3f} ± {row['std']:.3f}", axis = 1)
