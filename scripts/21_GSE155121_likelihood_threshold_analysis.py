##########################################################
## ccAFv2:  Likelihood threshold analysis GSE155121     ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccnn

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

_classifier = keras.models.load_model('ccAFv2_model.h5')
_genes = list(pd.read_csv('ccAFv2_genes.csv', index_col=0, header=0)['human_ensembl'])
_classes = list(pd.read_csv('ccAFv2_classes.txt', header=None)[0])

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
tags = ['GSE155121']
tag = 'NSC'
ws = ['W3-1', 'W4-1', 'W4-2', 'W4-3', 'W5-1', 'W5-2', 'W5-3', 'W6-1', 'W7-1', 'W8-1', 'W9-1', 'W9-2', 'W12-1']
resdir = 'data'
output = 'results/cutoff_analysis'
if not os.path.exists(output):
    os.makedirs(output)
    print ("Directory created")
else:
    print("Directory already exists")

# Common parameters
cutoffs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
nfolds = 10

# Load up data files and insert necessary information into objects before merge
datasets = {}
for tag1 in tags:
    for ws1 in ws:
        print('\nLoading '+tag1+' '+ws1+' data...')
        if tag1 == 'GSE155121':
            resdir2 = resdir+'/'+tag1+'/NSC'
            datasets[tag1] = sc.read_h5ad(resdir2+'/'+ws1+'_normalized_ensembl.h5ad')
            datasets[tag1].obs['new_clusters'] = datasets[tag1].obs['ccAFv2']
        datasets[tag1].obs['dataset'] = tag1
        datasets[tag1].obs['week_stage'] = ws1
        print(datasets[tag1].shape)
        #------------------------
        # Cross validation
        #-----------------------
        # Save as data as new variable
        ccAF1_scanpy = datasets[tag1]
        # Initialize helper vars/indices for subsetting genes
        nCells = ccAF1_scanpy.shape[0]
        allInds = np.arange(0, nCells)
        numSamples = round(0.9*nCells)
        # Make folder to store downstream results
        savedir = resdir2+'/'+ws1+'/'+save_fold
        if not os.path.exists(savedir):
            os.makedirs(savedir)
            print ("Directory created")
        else:
            print("Directory already exists")
        # Cross validation within dataset
        ccseurat_lab = ccAF1_scanpy.obs['Phase']
        for cutoff1 in cutoffs:
            print('\nUnknown cutoff...', cutoff1)
            barcodes = []
            pred = []
            truelab = []
            seuratlab = []
            for k in range(nfolds):
                samp1 = sample_without_replacement(nCells, numSamples, random_state = 1234 + k)
                print('Predicting...')
                testPredLbls = predict_labels(ccAF1_scanpy[samp1,:], cutoff = cutoff1)
                trueLbls = ccAF1_scanpy[samp1,:].obs['new_clusters']
                seuratLbls = ccAF1_scanpy[samp1,:].obs['Phase']
                barcodesLbls = ccAF1_scanpy[samp1,:].obs_names
                truelab.extend(trueLbls)
                seuratlab.extend(seuratLbls)
                pred.extend(testPredLbls[0])
                barcodes.extend(barcodesLbls)
            # Dataframe of true labels, predictions, probabilities for all iterations
            ind1 = pd.DataFrame(barcodes).rename(columns={0:'cell_barcodes'})
            col1 = pd.DataFrame(seuratlab).rename(columns={0:'ccSeurat'})
            col2 = pd.DataFrame(truelab).rename(columns={0:'True Labels'})
            col3 = pd.DataFrame(pred).rename(columns = {0:'Predictions'})
            DF = pd.concat([col1, col2, col3], axis=1)
            DF.index = ind1['cell_barcodes']
            sum(DF['True Labels'] == DF['Predictions'])/len(DF)
            # Save out data
            print('\nSaving out CV sensitivity analysis file...')
            DF.to_csv(savedir+'/'+ws1+'_ccAFv2_CV_sensitivity_analysis_'+str(cutoff1)+'.csv')

### All week stages together
# Adjusted Mutual Score
input = {}
results = {}
cells_predicted = {}
for tag1 in tags:
    input[tag1] = {}
    results[tag1] = {}
    cells_predicted[tag1] = {}
    for ws1 in ws:
        print('Setting file path: GSE155121 NSC '+ws1+'\n')
        resdir2 = resdir+'/'+tag1+'/NSC'
        savedir = resdir2+'/'+ws1+'/'+save_fold
        input[tag1][ws1] = {}
        results[tag1][ws1] = {}
        cells_predicted[tag1][ws1] = {}
        for cutoff1 in cutoffs:
            print('Reading in file:\nCV sensitivity analysis for unknown cutoff '+str(cutoff1)+'\n')
            input[tag1][ws1][cutoff1] = pd.read_csv(savedir+'/'+ws1+'_ccAFv2_CV_sensitivity_analysis_'+str(cutoff1)+'.csv', low_memory=False)
            # Set Unknowns to 'NaN' so can remove them easily during downstream analysis
            input[tag1][ws1][cutoff1]['Predictions'].replace('Unknown', np.nan, inplace=True)
        # we are using cutoff 0.0 but could easily use any cutoff
        numSamples = (len(input[tag1][ws1][0.0]))/nfolds
        for cutoff1 in cutoffs:
            results[tag1][ws1][cutoff1] = []
            cells_predicted[tag1][ws1][cutoff1] = []
            seuratlab = input[tag1][ws1][cutoff1]['ccSeurat']
            predictlab = input[tag1][ws1][cutoff1]['Predictions']
            for k in range(nfolds):
                bind = int(numSamples*k)
                eind = int(numSamples*(k+1))
                truelab = seuratlab.iloc[bind:eind]
                predlab = predictlab.iloc[bind:eind]
                if len(predlab.dropna()) <= 20:
                    results[tag1][ws1][cutoff1].append(np.nan)
                    cells_predicted[tag1][ws1][cutoff1].append(np.nan)
                else:
                    results[tag1][ws1][cutoff1].append(adjusted_mutual_info_score(truelab[predlab.dropna().index], predlab.dropna()))
                    # Cells predicted
                    cells_predicted[tag1][ws1][cutoff1].append(1 - sum(predlab.isna())/len(predlab))
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    for ws1 in ws:
            for cutoff1 in cutoffs:
                df = pd.merge(pd.DataFrame(range(nfolds), results[tag1][ws1][cutoff1]).reset_index().rename(columns = {'index':'adjusted_mutual_score'}), pd.DataFrame(range(nfolds), cells_predicted[tag1][ws1][cutoff1]).reset_index().rename(columns = {'index':'cells_predicted'}), on=0)
                df['week_stage'] = ws1
                df['cutoff'] = cutoff1
                col1.append(list(df['adjusted_mutual_score']))
                col2.append(list(df['cells_predicted']))
                col3.append(list(df[0]))
                col4.append(list(df['week_stage']))
                col5.append(list(df['cutoff']))
    col1 = [item for sublist in col1 for item in sublist]
    col2 = [item for sublist in col2 for item in sublist]
    col3 = [item for sublist in col3 for item in sublist]
    col4 = [item for sublist in col4 for item in sublist]
    col5 = [item for sublist in col5 for item in sublist]
    df2 = pd.DataFrame([col1, col2, col3, col4, col5]).T
    df2.rename(columns={0:'adjusted_mutual_score', 1:'cells_predicted', 2:'k', 3:'week_stage', 4:'cutoff'}, inplace=True)
    # Save out all information
    df2.to_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis.csv')
    all1 = df2
    ttest_sig = {}
    for cutoff1 in cutoffs[1:]:
        ref = all1[all1['cutoff']==0.0]['adjusted_mutual_score']
        comp = all1[all1['cutoff']==cutoff1]['adjusted_mutual_score']
        ttest_sig[cutoff1] = scipy.stats.ttest_ind(ref, comp)[1]
    tmp = pd.DataFrame.from_dict(ttest_sig, orient='index') # significant at 0.5!!!!!!
    tmp.to_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_statistics.csv')
    ttest_sig_corrected = statsmodels.stats.multitest.fdrcorrection(list(ttest_sig.values()))[1]
    tmp2 = pd.DataFrame(ttest_sig_corrected)
    tmp2.index = tmp.index
    tmp2.to_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_statistics_fdr_corrected.csv')
    for plot1 in ['adjusted_mutual_score', 'cells_predicted']:
        sns.set(style="whitegrid", font_scale=3)
        fig, ax1 = plt.subplots(figsize=(40,20))
        sns.boxplot(hue = "week_stage", y = plot1, x = "cutoff", data = all1, ax = ax1, palette = "husl")
        sns.move_legend(ax1, "best", bbox_to_anchor=(1,1), ncol=1)
        ax1.set(ylabel=plot1)
        plt.savefig(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_'+plot1+'.pdf')
        plt.clf()
        # Plot median scatter plots
        median1 = {}
        for ws1 in ws:
            median1[ws1] = []
            for cutoff1 in cutoffs:
                median1[ws1].append(np.median(all1[all1['week_stage'] == ws1][all1[all1['week_stage'] == ws1]['cutoff']==cutoff1][plot1]))
        all_ws_medians = pd.DataFrame(median1)
        all_ws_medians.index = cutoffs
        all_ws_medians.index.name = 'cutoff'
        all_ws_medians.to_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_median_'+plot1+'_cutoff_for_all_week_stages.csv')
        plt.clf()

# Plotting
amis = pd.read_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_median_adjusted_mutual_score_cutoff_for_all_week_stages.csv', index_col=0)
cps = pd.read_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_median_cells_predicted_cutoff_for_all_week_stages.csv', index_col=0)
median_amis = []
median_cps = []
for cutoff1 in cutoffs:
    median_amis.append(np.median(amis.loc[cutoff1]))
    median_cps.append(np.median(cps.loc[cutoff1]))

x = np.array(median_cps)*100
y = np.array(median_amis)
descrip = np.array(amis.index)
sns.set(style="whitegrid", font_scale=0.5)
fig, ax = plt.subplots()
ax.errorbar(x, y, fmt='o', ecolor = 'black')
ax.set_xlabel('Cells Predicted (%)')
ax.set_ylabel('Adjusted Mutual Information Score')
ax.set_ylim(0.3,0.5)
ax.set_xlim(50,100)
for i, txt in enumerate(descrip):
    ax.annotate(txt, (x[i], y[i]))
plt.savefig(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_AMI_cellsPredicted_together_all_week_stages_together.pdf')
plt.clf()
