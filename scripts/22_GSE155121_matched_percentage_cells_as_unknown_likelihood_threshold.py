##########################################################
## ccAFv2:  Matched unknown cell loss GSE155121         ##
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

#--------------------------------
# Set up section / load packages
#--------------------------------
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

# Folder set up
tags = ['GSE155121']
tag = 'NSC'
ws = ['W3-1', 'W4-1', 'W4-2', 'W4-3', 'W5-1', 'W5-2', 'W5-3', 'W6-1', 'W7-1', 'W8-1', 'W9-1', 'W9-2', 'W12-1']
resdir = 'data'
output = 'results/cutoff_analysis'
savedir = 'results/matched_cell_perc_analysis'
if not os.path.exists(savedir):
    os.makedirs(savedir)



# Load cutoff analysis info
cells_pred_cutoff = pd.read_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_median_cells_predicted_cutoff_for_all_week_stages.csv', index_col=0)
adjusted_mut_cutoff = pd.read_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis_median_adjusted_mutual_score_cutoff_for_all_week_stages.csv', index_col=0)
all_info_cutoff = pd.read_csv(output+'/GSE155121_NSCs_ccAFv2_likelihood_threshold_analysis.csv', index_col=0)

#----------------
# Load data
#----------------
nfolds = 10
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
        #------------------------
        # Cross validation
        #-----------------------
        # Save as data as new variable
        ccAF1_scanpy = datasets[tag1]
        # Initialize helper vars/indices for subsetting genes
        nCells = ccAF1_scanpy.shape[0]
        # Cross validation within dataset
        ccseurat_labs = ccAF1_scanpy.obs['Phase']
        truelabs = ccAF1_scanpy.obs['new_clusters']
        # Random chance?
        cell_percs = round(cells_pred_cutoff[ws1],2).dropna()
        chance_pred = dict(zip(cell_percs,[[] for i in range(len(cell_percs))]))
        ccseurat_lab = dict(zip(cell_percs,[[] for i in range(len(cell_percs))]))
        truelab = dict(zip(cell_percs,[[] for i in range(len(cell_percs))]))
        for cell_perc1 in [i for i in list(chance_pred.keys()) if i != 0]:
            print('\n' +str(cell_perc1))
            for k in range(nfolds):
                samp2 = sample_without_replacement(nCells, int(nCells*cell_perc1), random_state = 1234 + k)
                samp3 = pd.Series(ccAF1_scanpy.obs_names).iloc[samp2]
                ccseurat_lab[cell_perc1].extend(ccAF1_scanpy[list(samp3),:].obs['Phase'])
                truelab[cell_perc1].extend(ccAF1_scanpy[list(samp3),:].obs['new_clusters'])
                print('Predicting...')
                testPredLbls = predict_labels(ccAF1_scanpy[list(samp3),:], cutoff=0.0)
                chance_pred[cell_perc1].extend(testPredLbls[0])
            # Dataframe of true labels, predictions, probabilities for all iterations
            col1 = pd.DataFrame(ccseurat_lab[cell_perc1]).rename(columns={0:'ccSeurat'})
            col2 = pd.DataFrame(truelab[cell_perc1]).rename(columns={0:'True Labels'})
            col3 = pd.DataFrame(chance_pred[cell_perc1]).rename(columns={0:'Pred Labels'})
            DF = pd.concat([col1, col2, col3], axis=1)
            # Save out data
            print('\nSaving out CV sensitivity analysis file...')
            DF.to_csv(savedir+'/'+ws1+'_ccAFv2_CV_matched_'+str(cell_perc1)+'_perc_analysis.csv')

#------------------------
# Downstream analysis
#------------------------
input = {}
results = {}
cells_predicted = {}
for tag1 in tags:
    input[tag1] = {}
    results[tag1] = {}
    cells_predicted[tag1] = {}
    for ws1 in ws:
        print('Setting file path: GSE155121 NSC '+ws1+'\n')
        input[tag1][ws1] = {}
        results[tag1][ws1] = {}
        cells_predicted[tag1][ws1] = {}
        cell_percs = round(cells_pred_cutoff[ws1],2).dropna()
        cell_percs_run = [i for i in cell_percs.unique() if i != 0]
        for cell_perc1 in cell_percs_run:
            print('\n' +str(cell_perc1))
            results[tag1][ws1][cell_perc1] = []
            cells_predicted[tag1][ws1][cell_perc1] = []
            input[tag1][ws1][cell_perc1] = pd.read_csv(savedir+'/'+ws1+'_ccAFv2_CV_matched_'+str(cell_perc1)+'_perc_analysis.csv', low_memory=False)
            numSamples = len(input[tag1][ws1][cell_perc1])/nfolds
            truelabs = input[tag1][ws1][cell_perc1]['ccSeurat']
            # Separate by nfolds so can get variance
            for k in range(nfolds):
                bind = int(numSamples*k)
                eind = int(numSamples*(k+1))
                truelab = truelabs.iloc[bind:eind]
                predlab = input[tag1][ws1][cell_perc1].iloc[bind:eind]['Pred Labels']
                if len(predlab.dropna()) <= 20:
                    results[tag1][ws1][cell_perc1].append(np.nan)
                    cells_predicted[tag1][ws1][cell_perc1].append(np.nan)
                else:
                    # Adjusted mutual score; drop unknowns
                    results[tag1][ws1][cell_perc1].append(adjusted_mutual_info_score(truelab[predlab.dropna().index], predlab.dropna()))
                    # Cells predicted
                    cells_predicted[tag1][ws1][cell_perc1].append(1 - sum(predlab.isna())/len(predlab))
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    for ws1 in ws:
        cell_percs = round(cells_pred_cutoff[ws1],2).dropna()
        cell_percs_run = [i for i in cell_percs.unique() if i != 0]
        for cell_perc1 in cell_percs_run:
            df = pd.merge(pd.DataFrame(range(nfolds), results[tag1][ws1][cell_perc1]).reset_index().rename(columns = {'index':'adjusted_mutual_score'}), pd.DataFrame(range(nfolds), cells_predicted[tag1][ws1][cell_perc1]).reset_index().rename(columns = {'index':'cells_predicted'}), on=0)
            df['week_stage'] = ws1
            df['percent_cells'] = cell_perc1
            col1.append(list(df['adjusted_mutual_score']))
            col2.append(list(df['cells_predicted']))
            col3.append(list(df[0]))
            col4.append(list(df['percent_cells']))
            col5.append(list(df['week_stage']))
    col1 = [item for sublist in col1 for item in sublist]
    col2 = [item for sublist in col2 for item in sublist]
    col3 = [item for sublist in col3 for item in sublist]
    col4 = [item for sublist in col4 for item in sublist]
    col5 = [item for sublist in col5 for item in sublist]
    df2 = pd.DataFrame([col1, col2, col3, col4, col5]).T
    df2.rename(columns={0:'adjusted_mutual_score', 1:'cells_predicted', 2:'k', 3:'percent_cells', 4:'week_stage'}, inplace=True)
    df2.to_csv(savedir+'/GSE155121_NSCs_matched_percent_cells_for_each_week_stage.csv')
    median1 = {}
    median_cutoff = {}
    ttest = {}
    matched_cutoff_ttest = {}
    compare_box_plot_comp = {}
    compare_box_plot_ref = {}
    for ws1 in ws:
        # Find medians for both cutoff analysis and random cell loss analysis
        median1[ws1] = {}
        median_cutoff[ws1] = {}
        ttest[ws1] = {}
        matched_cutoff_ttest[ws1] = {}
        compare_box_plot_comp[ws1] = {}
        compare_box_plot_ref[ws1] = {}
        cell_percs = round(cells_pred_cutoff[ws1],2).dropna()
        cell_percs_run = [i for i in cell_percs.unique() if i != 0]
        for cell_perc1 in cell_percs_run:
            comp = df2[df2['week_stage']==ws1][df2[df2['week_stage'] == ws1]['percent_cells']==cell_perc1]['adjusted_mutual_score']
            find_cutoff = pd.DataFrame(round(cells_pred_cutoff[ws1], 2)).loc[pd.DataFrame(round(cells_pred_cutoff[ws1], 2))[ws1] == cell_perc1]
            cutoff_val = list(find_cutoff.index).pop()
            ws_cutoff = all_info_cutoff[all_info_cutoff['week_stage']==ws1][all_info_cutoff[all_info_cutoff['week_stage']==ws1]['cutoff'] == cutoff_val]
            median_cutoff[ws1][cell_perc1] = np.median(ws_cutoff['adjusted_mutual_score'])
            ref = ws_cutoff['adjusted_mutual_score']
            ttest[ws1][cell_perc1] = scipy.stats.ttest_ind(ref, comp)[1]
            matched_cutoff_ttest[ws1][cutoff_val] = scipy.stats.ttest_ind(ref, comp)[1]
            compare_box_plot_comp[ws1][cutoff_val] = list(comp.values)
            compare_box_plot_ref[ws1][cutoff_val] = list(ref.values)
    # Save out ttest statistics
    df3 = pd.DataFrame(matched_cutoff_ttest)
    df3.to_csv(savedir+'/GSE155121_NSC_adjusted_mutual_info_score_ttest_values_between_unknown_cutoff_and_matched_random_cells_dropped.csv')
    df3_corrected = {}
    for ws1 in ws:
        df3_corrected[ws1] = statsmodels.stats.multitest.fdrcorrection(df3[ws1])[1]
    df4 = pd.DataFrame(df3_corrected)
    df4.index = df3.index
    df4.to_csv(savedir+'/GSE155121_NSC_adjusted_mutual_info_score_fdr_corrected_ttest_values_between_unknown_cutoff_and_matched_random_cells_dropped.csv')
    # Plot all together in one PDF
    with PdfPages(r'results/matched_cell_perc_analysis/GSE155121_NSC_dropout_matched_perc_cells_all_together.pdf') as export_pdf:
        for ws1 in ws:
            tmp1 = pd.DataFrame(compare_box_plot_comp[ws1]).melt()
            tmp1['week_stage'] = ws1
            tmp1['analysis'] = 'matched_percentage_cell_loss'
            #all_comp = all_comp.append(tmp1)
            tmp2 = pd.DataFrame(compare_box_plot_ref[ws1]).melt()
            tmp2['week_stage'] = ws1
            tmp2['analysis'] = 'unknown_cutoff'
            to_plot = tmp1.append(tmp2)
            sns.set(style="whitegrid", font_scale=3)
            fig, ax1 = plt.subplots(figsize=(40,20))
            sns.boxplot(hue = 'analysis', y = "value", x = "variable", data = to_plot, ax = ax1, palette = "husl")
            ax1.set(ylabel="adjusted_mutual_info_score", xlabel = 'cutoff')
            ax1.set_title('GSE155121_'+ws1)
            #ax1.set_ylim(0,1)
            plt.savefig(savedir+'/GSE155121_likelihood_threshold_versus_random_cells_dropped_'+ws1+'.pdf')
            export_pdf.savefig(fig)
