##########################################################
## ccAFv2:  testing ccAFv2 R vs. python                 ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

#docker run -it -v '/home/soconnor/old_home:/files' cplaisier/ccnn

# General
from importlib.resources import path
import numpy as np
import pandas as pd
import os
from scipy.sparse import isspmatrix
import scanpy as sc
sc.settings.verbosity = 0
import calendar
import time
time_stamp = calendar.timegm(time.gmtime())

# ccAFv2
from sklearn.preprocessing import StandardScaler, LabelEncoder
import tensorflow as tf
from tensorflow import keras

# Stop warning messages for cudart
import logging
logging.disable(logging.WARNING)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "3"
logging.getLogger('tensorflow').disabled = True

################
## Load model ##
################
os.chdir('ccAFv2_py')

with path('ccAF', 'ccAFv2_model.h5') as inPath:
    _classifier = keras.models.load_model(inPath)
with path('ccAF', 'ccAFv2_genes.csv') as inPath:
    _genes = list(pd.read_csv(inPath, index_col=0, header=0)['human_ensembl'])
with path('ccAF', 'ccAFv2_classes.txt') as inPath:
    _classes = list(pd.read_csv(inPath, header=None)[0])

###############
## Functions ##
###############

def _scale(data):
    """
    Standardize or normalize numeric data using scikit-learn's StandardScaler.

    Parameters:
    - data: 2D NumPy array or list of lists

    Returns:
    - Scaled data (NumPy array)
    """
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    return scaled_data

# Prepare test data for predicting
def _prep_predict_data(data, genes):
    """
    prep_predict_data takes in a pandas dataframe and the trained ccAFv2 model.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame of scRNA-seq data to be classified.
    model : keras.models.sequential
        Trained ccAFv2 sequential keras model.

    Returns
    -------
    pd.Series
        Series of labels for each single cell.

    """
    # Remove all genes with zero counts
    data.var_names_make_unique()
    sc.pp.filter_genes(data, min_cells=1)
    # Restrict to classifier genes
    data2 = data[:,list(set(genes).intersection(data.var_names))]
    # Scale data
    if isspmatrix(data.X):
        data2 = pd.DataFrame(data2.X.todense(), index = data2.obs_names, columns = data2.var_names)
    else:
        data2 = pd.DataFrame(data2.X, index = data2.obs_names, columns = data2.var_names)
    data3 = pd.DataFrame(_scale(data2), index = data2.index, columns = data2.columns)
    # Add minimum values for missing genes
    missing = set(genes).difference(data3.columns)
    if len(missing)>0:
        data4 = pd.concat([data3, pd.DataFrame(data3.values.min(), index=data3.index, columns = missing)], axis=1)
        return data4[list(genes)]
    else:
        return data3

# Predict labels with rejection
def predict_labels(new_data, classifier=_classifier, genes=_genes, classes=_classes, cutoff=0.5):
    """
    predict_new_data takes in a pandas dataframe and the trained ccAFv2 model.

    Parameters
    ----------
    new_data : annData object
         of scRNA-seq data to be classified.
    cutoff : float
        The cutoff for likelihoods from the neural network classifier model.

    Returns
    -------
    pd.Series
        Series of labels for each single cell.

    """
    pred_data = _prep_predict_data(new_data, genes)
    probabilities = _predict_new_data(pred_data, classifier)
    labels = np.array([classes[np.argmax(i)] for i in probabilities])
    labels[np.where([np.max(i) < cutoff for i in probabilities])] = 'Unknown'
    return labels, probabilities

# Predict ccAFv2 labels for new data
def _predict_new_data(new_data, classifier):
    """
    predict_new_data takes in a pandas dataframe and the trained ccAFv2 model.

    Parameters
    ----------
    new_data : pd.DataFrame
        DataFrame of scRNA-seq data to be classified.
    model : keras.models.sequential
        Trained ccAFv2 sequential keras model.

    Returns
    -------
    pd.Series
        Series of labels for each single cell.

    """
    return classifier.predict(new_data)

# Folder and tag set up
# Test GSE155121 W8-1
os.chdir('../ccNN')
tags = ['GSE155121']
tag = 'NSC'
ws = 'W8-1'
resdir = 'testData'
resdir5 = 'results'

# Read in data
data1 = sc.read_h5ad(resdir+'/GSE155121/NSC/'+ws+'_normalized_ensembl_test2.h5ad') # has ccAFv2 R in metadata
predictions = predict_labels(data1)
data1.obs['ccAFv2_py'] = predictions[0]

# Do they match?
sum(data1.obs['ccAFv2'] == data1.obs['ccAFv2_py'])/len(data1) # 1.0  yes

# Save out obs to csv
pd.DataFrame(data1.obs).to_csv('R_vs_python_GSE155121_compare.csv')
