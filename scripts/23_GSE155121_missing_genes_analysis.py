##########################################################
## ccAFv2: GSE155121 sensitivityAnalysis.py             ##
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

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim
from sklearn.metrics import adjusted_mutual_info_score

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
import scipy
from scipy import stats
import math
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

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
# change directory for testing
os.chdir('../ccNN')
tags = ['GSE155121']
tag = 'NSC'
ws = ['W8-1']
resdir = 'testData'
resdir5 = 'results'

# Common parameters
cutoffs = [0.5]
nfolds = 10

# Load up data files and insert necessary information into objects before merge
datasets = {}
for tag1 in tags:
    for ws1 in ws:
        print('\nLoading '+tag1+' '+ws1+' data...')
        if tag1 == 'GSE155121':
            resdir7 = resdir+'/'+tag1+'/NSC'
            datasets[tag1] = sc.read_h5ad(resdir+'/GSE155121/NSC/'+ws1+'_normalized_ensembl_test2.h5ad')
        datasets[tag1].obs['dataset'] = tag1
        datasets[tag1].obs['week_stage'] = ws1
        print(datasets[tag1].shape)

#------------------------
# Cross validation
#-----------------------

# Save as data as new variable
ccAF1_scanpy = datasets[tag1]
ccAF1_scanpy.obs['new_clusters'] = ccAF1_scanpy.obs['ccAFv2']
# Parameters for CV
nfolds = 10
# Initialize helper vars/indices for subsetting genes
nCells = ccAF1_scanpy.shape[0]
nGenes = ccAF1_scanpy.shape[1]
GeneInds = np.arange(0, nGenes)

#################
### ccAFv2 CV ###
#################

# Make folder to store downstream results
if not os.path.exists('results/SensitivityAnalysis'):
    os.makedirs('results/SensitivityAnalysis')
    print ("Directory created")
else:
    print("Directory already exists")

if not os.path.exists('results/SensitivityAnalysis/final_cutoff'):
    os.makedirs('results/SensitivityAnalysis/final_cutoff')
    print ("Directory created")
else:
    print("Directory already exists")


# Cross validation within dataset
# Percentage of genes to test for cross-validation classification
cutoffs = [0.5, 0.7, 0.9]
percs = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
for cutoff1 in cutoffs:
    print('\nUnknown cutoff...', cutoff1)
    truelab = ccAF1_scanpy.obs['new_clusters']
    pred = dict(zip(percs,[[] for i in range(len(percs))]))
    print('\nccAFv2 cross-validation...')
    for perc1 in percs:
        print('\nPercentage of genes in dataset:', perc1)
        for k in range(nfolds):
            print('ccAFv2 round '+str(k)+'...')
            samp2 = sample_without_replacement(len(_genes), int(len(_genes)*perc1), random_state = 1234 + k)
            samp3 = pd.Series(_genes).iloc[samp2]
            #truelab.extend(ccAF1_scanpy.obs['new_clusters'])
            print('Predicting...')
            testPredLbls = predict_labels(ccAF1_scanpy[:,ccAF1_scanpy.var_names.intersection(list(samp3))], cutoff = cutoff1)
            pred[perc1].extend(testPredLbls[0])
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.concat([pd.DataFrame(list(truelab)*10).rename(columns={0:'True Labels'}), pd.DataFrame(pred)], axis=1)[['True Labels']+percs]
    #DF.replace('nan', np.nan, inplace=True)
    # Save out data
    DF.to_csv('results/SensitivityAnalysis/final_cutoff/'+tag1+'_ccAFv2_CV_sensitivity_analysis_'+str(cutoff1)+'_032224.csv')


# Read in data for plotting
resdir9 = 'results/SensitivityAnalysis/final_cutoff'
nfolds = 10
labs = ccAF1_scanpy.obs['Phase']
extra_col = list(labs)*10
#cutoffs = [0.5, 0.7, 0.9]
cutoffs = [0.9]
percs = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
#percs = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2] #0.9
input = {}
results = {}
cells_predicted = {}
for cutoff1 in cutoffs:
    file1 = pd.read_csv(resdir9+'/'+tag1+'_ccAFv2_CV_sensitivity_analysis_'+str(cutoff1)+'_032224.csv')
    file1['ccSeurat'] = extra_col
    input[cutoff1] = file1
    numSamples = len(input[cutoff1])/nfolds
    results[cutoff1] = {}
    cells_predicted[cutoff1] = {}
    for perc1 in percs:
        results[cutoff1][perc1] = []
        defined_cell_states = input[cutoff1]['True Labels']
        truelabs = input[cutoff1]['ccSeurat']
        cells_predicted[cutoff1][perc1] = []
        for k in range(nfolds):
            bind = int(numSamples*k)
            eind = int(numSamples*(k+1))
            defined_cell_state = defined_cell_states.iloc[bind:eind]
            truelab = truelabs.iloc[bind:eind] # seurat labels
            #truelab = defined_cell_states.iloc[bind:eind] # ccAF labels
            predlab = input[cutoff1][str(perc1)].iloc[bind:eind]
            # Change 'Unknown' to NaN so can drop those out
            predlab.replace('Unknown', np.nan, inplace=True)
            # Adjusted mutual score; drop unknowns
            results[cutoff1][perc1].append(adjusted_mutual_info_score(truelab[predlab.dropna().index], predlab.dropna()))
            # Cells predicted
            cells_predicted[cutoff1][perc1].append(1 - sum(predlab.isna())/len(predlab))
col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
for cutoff1 in cutoffs:
    for perc1 in percs:
        df = pd.merge(pd.DataFrame(range(nfolds), results[cutoff1][perc1]).reset_index().rename(columns = {'index':'adjusted_mutual_score'}), pd.DataFrame(range(nfolds), cells_predicted[cutoff1][perc1]).reset_index().rename(columns = {'index':'cells_predicted'}), on=0)
        df['cutoff'] = cutoff1
        df['percentage_genes'] = perc1
        col1.append(list(df['adjusted_mutual_score']))
        col2.append(list(df['cells_predicted']))
        col3.append(list(df[0]))
        col4.append(list(df['cutoff']))
        col5.append(list(df['percentage_genes']))
col1 = [item for sublist in col1 for item in sublist]
col2 = [item for sublist in col2 for item in sublist]
col3 = [item for sublist in col3 for item in sublist]
col4 = [item for sublist in col4 for item in sublist]
col5 = [item for sublist in col5 for item in sublist]
df2 = pd.DataFrame([col1, col2, col3, col4, col5]).T
df2.rename(columns={0:'adjusted_mutual_score', 1:'cells_predicted', 2: 'k', 3: 'cutoff', 4:'percentage_genes'}, inplace=True)
df2.to_csv(resdir9+'/'+tag1+'_percentage_genes_remove_unknowns_'+str(cutoff1)+'_032224.csv')
# Together plots
sns.set(style="whitegrid", font_scale=3)
fig, ax1 = plt.subplots(figsize=(40,20))
df3 = df2.melt(id_vars = ['percentage_genes', 'k', 'cutoff'], var_name='metric', value_name='value')
sns.boxplot(hue = "percentage_genes", y = "value", x = "metric", data = df3, ax = ax1, palette = "husl")
sns.move_legend(ax1, "best", bbox_to_anchor=(1,1), ncol=1)
ax1.set(ylabel="value")
ax1.set_ylim(0,1)
plt.savefig(resdir9+'/'+tag1+'_percentage_genes_boxplot_together_remove_unknowns_'+str(cutoff1)+'_032224.pdf')
plt.clf()


# Find metric medians
median_cp = {}
median_ami = {}
for perc1 in percs:
    median_cp[perc1] = np.median(df2[df2['percentage_genes'] == perc1]['cells_predicted'])*100
    median_ami[perc1] = np.median(df2[df2['percentage_genes'] == perc1]['adjusted_mutual_score'])


df3 = pd.DataFrame(
    {'cells_predicted': list(median_cp.values()),
     'ami': list(median_ami.values())
    })
df3.index = percs
df3.to_csv(resdir9+'/'+tag1+'_all_metrics_percentage_genes_sensitivity_medians_'+str(cutoff1)+'_032224.csv')
