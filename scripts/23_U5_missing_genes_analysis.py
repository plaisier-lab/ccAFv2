##########################################################
## ccAFv2: U5 sensitivityAnalysis.py                    ##
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

def _scale(df1):
    """
    scale takes in a pandas dataframe and applies scales the values into Z-scores across rows.

    Parameters
    ----------
    df1 : pd.DataFrame
        DataFrame of scRNA-seq data to be scaled.

    Returns
    -------
    pd.DataFrame
        DataFrame of scRNA-seq data that has been scaled.
    """
    return (df1.subtract(df1.mean(axis=1),axis=0)).div(df1.std(axis=1),axis=0)

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
sets = ['U5']
p_cutoffs = {'U5':0.05}
resdir = 'data/'+sets[0]
savedir = 'results/missing_genes_analysis'
if not os.path.exists(savedir):
    os.makedirs(savedir)

# Load up conversion file
symEnsmbl = pd.read_csv('geneConversion/U5/filtered_feature_bc_matrix/features.tsv.gz', header = None, index_col=1, sep='\t')
tmp1 = pd.Series(symEnsmbl.index)
tmp1.loc[symEnsmbl.index.duplicated()] = [i+'.1' for i in symEnsmbl.loc[symEnsmbl.index.duplicated()].index]
symEnsmbl.index = pd.Index(tmp1)

# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv('geneConversion/Whitfield/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
hgncEnsmbl = hgncEnsmbl.loc[~hgncEnsmbl['Ensembl ID(supplied by Ensembl)'].isnull()]

ensmblHgnc = pd.Series(hgncEnsmbl.index)
ensmblHgnc.index = list(hgncEnsmbl['Ensembl ID(supplied by Ensembl)'])

hgncPrevEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Previous symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Previous symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncPrevEnsmbl[j] = ensmbl

hgncAliasEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Alias symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Alias symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncAliasEnsmbl[j] = ensmbl

#---------------------------------------------
# Load data and subset genes to top marker genes
#--------------------------------------------
datasets = {}
for set1 in sets:
    # Load in normalized U5-hNSC data
    print('\nLoading U5-hNSC scRNA-seq data...')
    datasets[set1] = sc.read_loom(resdir+'/'+set1+'_normalized_gene_symbols.loom')
    datasets[set1].obs['new_clusters'] = datasets[set1].obs['ccAF']
    datasets[set1].obs['dataset'] = set1
    # Convert genes to Ensembl
    print('\nConverting gene symbols to Ensembl...')
    common = [i for i in datasets[set1].var_names if i in symEnsmbl[0]]
    datasets[set1]= datasets[set1][:,common]
    datasets[set1].var_names = pd.Index(symEnsmbl.loc[datasets[set1].var_names,0], name='Ensembl')
    # Subset marker genes
    print('\nSubsetting data genes to marker genes...')
    mgenes = pd.read_csv('markerGenes/U5_scTransform_Markers_together.csv', header = 0, index_col = 0)
    # Filter marker genes by avg_log2FC
    mgenes = mgenes.loc[mgenes['avg_log2FC'].ge(0.25),:]
    # Remove G1/other cells for ccAFv2 training
    mgenes = mgenes[mgenes['cluster'] != 'G1/other']
    # Filter marker genes by p_val_adj
    mgenes = mgenes.loc[mgenes['p_val_adj'].le(p_cutoffs[set1]),:]
    mgenes = list(mgenes['gene'])
    mgenes = list(set(mgenes))
    tmp1 = []
    missed = []
    for j in mgenes:
        if j in hgncEnsmbl.index:
            tmp1.append(hgncEnsmbl.loc[j,'Ensembl ID(supplied by Ensembl)'])
        elif j in hgncPrevEnsmbl:
            tmp1.append(hgncPrevEnsmbl[j])
        elif j in hgncAliasEnsmbl:
            tmp1.append(hgncAliasEnsmbl[j])
        else:
            missed.append(j)
    mg1 = tmp1
    mg1 = list(set(mg1).intersection(datasets[set1].var_names))
    print('U5 marker genes: '+str(len(mg1)))
    # Subset data to marker genes
    datasets[set1] = datasets[set1][:,mg1]

#------------------------
# Cross validation
#-----------------------

# Save as data as new variable
ccAF1_scanpy = datasets[set1]
ccAFv2_calls = pd.read_csv(resdir+'/'+set1+'_ccAFv2_calls.csv', index_col = 0)
# Set column name
ccAF1_scanpy.obs['ccAFv2'] = ccAFv2_calls['x']
ccAF1_scanpy.obs['new_clusters'] = ccAF1_scanpy.obs['ccAF']
# Parameters for CV
nfolds = 10
# Initialize helper vars/indices for subsetting genes
nCells = ccAF1_scanpy.shape[0]
nGenes = ccAF1_scanpy.shape[1]
GeneInds = np.arange(0, nGenes)

#################
### ccAFv2 CV ###
#################

# Cross validation within dataset
# Percentage of genes to test for cross-validation classification
cutoffs = [0.5]
#cutoffs = [0.5, 0.7, 0.9]
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
            testPredLbls = predict_labels(ccAF1_scanpy[:,list(samp3)], cutoff = cutoff1)
            pred[perc1].extend(testPredLbls[0])
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.concat([pd.DataFrame(list(truelab)*10).rename(columns={0:'True Labels'}), pd.DataFrame(pred)], axis=1)[['True Labels']+percs]
    #DF.replace('nan', np.nan, inplace=True)
    # Save out data
    DF.to_csv(savedir+'/'+set1+'_ccAFv2_CV_missing_genes_analysis_'+str(cutoff1)+'.csv')


# Read in data for plotting
ccSeurat_calls = pd.read_csv(resdir+'/'+set1+'_ccSeurat_calls.csv')
ccSeurat_calls['x'].replace({'G2M': 'G2/M'}, inplace=True)
extra_col = list(ccSeurat_calls['x'])*10
cutoffs = [0.5]
percs = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
input = {}
results = {}
f1 = {}
cells_predicted = {}
error_rate = {}
for cutoff1 in cutoffs:
    file1 = pd.read_csv(savedir+'/'+set1+'_ccAFv2_CV_missing_genes_'+str(cutoff1)+'.csv')
    file1['ccSeurat'] = extra_col
    input[cutoff1] = file1
    numSamples = len(input[cutoff1])/nfolds
    results[cutoff1] = {}
    f1[cutoff1] = {}
    cells_predicted[cutoff1] = {}
    error_rate[cutoff1] = {}
    for perc1 in percs:
        results[cutoff1][perc1] = []
        defined_cell_states = input[cutoff1]['True Labels']
        truelabs = input[cutoff1]['ccSeurat']
        cells_predicted[cutoff1][perc1] = []
        error_rate[cutoff1][perc1] = []
        f1[cutoff1][perc1] = []
        for k in range(nfolds):
            bind = int(numSamples*k)
            eind = int(numSamples*(k+1))
            defined_cell_state = defined_cell_states.iloc[bind:eind]
            #truelab = truelabs.iloc[bind:eind] # seurat labels
            truelab = defined_cell_states.iloc[bind:eind] # ccAF labels
            predlab = input[cutoff1][str(perc1)].iloc[bind:eind]
            # Change 'Unknown' to NaN so can drop those out
            predlab.replace('Unknown', np.nan, inplace=True)
            # Adjusted mutual score; drop unknowns
            results[cutoff1][perc1].append(adjusted_mutual_info_score(truelab[predlab.dropna().index], predlab.dropna()))
            # Cells predicted
            cells_predicted[cutoff1][perc1].append(1 - sum(predlab.isna())/len(predlab))
            # Error Rate
            error_rate[cutoff1][perc1].append(1-sum(defined_cell_state[predlab.dropna().index] == predlab.dropna())/len(defined_cell_state[predlab.dropna().index])) if len(defined_cell_state[predlab.dropna().index]) != 0 else error_rate[cutoff1][perc1].append(0)
            # F1 score
            f1[cutoff1][perc1].append(classification_report(truelab[predlab.dropna().index], predlab.dropna(), output_dict=True, zero_division=0))
col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
for cutoff1 in cutoffs:
    for perc1 in percs:
        df = pd.merge(pd.merge(pd.DataFrame(range(nfolds), results[cutoff1][perc1]).reset_index().rename(columns = {'index':'adjusted_mutual_score'}), pd.DataFrame(range(nfolds), cells_predicted[cutoff1][perc1]).reset_index().rename(columns = {'index':'cells_predicted'}), on=0), pd.DataFrame(range(nfolds), error_rate[cutoff1][perc1]).reset_index().rename(columns = {'index':'error_rate'}), on=0)
        df['cutoff'] = cutoff1
        df['percentage_genes'] = perc1
        col1.append(list(df['adjusted_mutual_score']))
        col2.append(list(df['cells_predicted']))
        col3.append(list(df['error_rate']))
        col4.append(list(df[0]))
        col5.append(list(df['cutoff']))
        col6.append(list(df['percentage_genes']))
col1 = [item for sublist in col1 for item in sublist]
col2 = [item for sublist in col2 for item in sublist]
col3 = [item for sublist in col3 for item in sublist]
col4 = [item for sublist in col4 for item in sublist]
col5 = [item for sublist in col5 for item in sublist]
col6 = [item for sublist in col6 for item in sublist]
df2 = pd.DataFrame([col1, col2, col3, col4, col5, col6]).T
df2.rename(columns={0:'adjusted_mutual_score', 1:'cells_predicted', 2:'error_rate', 3: 'k', 4: 'cutoff', 5:'percentage_genes'}, inplace=True)
df2.to_csv(savedir+'/'+set1+'_percentage_genes_remove_unknowns_'+str(cutoff1)+'.csv')
# Together plots
sns.set(style="whitegrid", font_scale=3)
fig, ax1 = plt.subplots(figsize=(40,20))
df3 = df2.melt(id_vars = ['percentage_genes', 'k', 'cutoff'], var_name='metric', value_name='value')
sns.boxplot(hue = "percentage_genes", y = "value", x = "metric", data = df3, ax = ax1, palette = "husl")
sns.move_legend(ax1, "best", bbox_to_anchor=(1,1), ncol=1)
ax1.set(ylabel="value")
plt.savefig(savedir+'/'+set1+'_percentage_genes_boxplot_together_remove_unknowns_'+str(cutoff1)+'.pdf')
plt.clf()


# Find metric medians
median_er = {}
median_cp = {}
median_ami = {}
for perc1 in percs:
    median_er[perc1] = np.median(df2[df2['percentage_genes'] == perc1]['error_rate'])
    median_cp[perc1] = np.median(df2[df2['percentage_genes'] == perc1]['cells_predicted'])*100
    median_ami[perc1] = np.median(df2[df2['percentage_genes'] == perc1]['adjusted_mutual_score'])


df3 = pd.DataFrame(
    {'error_rate': list(median_er.values()),
     'cells_predicted': list(median_cp.values()),
     'ami': list(median_ami.values())
    })
df3.index = percs
df3.to_csv(savedir+'/U5_all_metrics_percentage_genes_sensitivity_medians_'+str(cutoff1)+'.csv')
