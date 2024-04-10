##########################################################
## ccAFv2:  CV_classification_methods_10_perc_holdout   ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccnn
#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/scrna_seq_velocity # ACTINN
#------------------------
# Set up / imports
#-----------------------

# General
import os
from os.path import exists
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

# Custom classes for classification
#import classifiersV3 as cl
import classifiersV6clp as cl

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap


###############
# CV analysis #
##############

# Folder set up
tags = ['U5']
cutoffs = {'U5':0.25}
resdir = 'data'
newfold = 'normalized'
resdir2 = 'geneConversion'
resdir5 = 'results'

# Load up data and insert necessary information into objects before merge
datasets = {}
for tag1 in tags:
    print('\nLoading '+tag1+' data...')
    if tag1 == 'U5':
        datasets[tag1] = sc.read_loom(resdir+'/'+newfold+'/final/U5_normalized_gene_symbols.loom')
        datasets[tag1].obs['new_clusters'] = datasets[tag1].obs['ccAF']
    if tag1 == 'ScienCell_hBMMSC':
        datasets[tag1] = sc.read_loom(resdir+'/'+newfold+'/'+tag1+'_normalized_new_labels_030423.loom')
    datasets[tag1].obs['dataset'] = tag1

# Load up conversion file
symEnsmbl = pd.read_csv(resdir2+'/'+tags[0]+'/filtered_feature_bc_matrix/features.tsv.gz', header = None, index_col=1, sep='\t')
tmp1 = pd.Series(symEnsmbl.index)
tmp1.loc[symEnsmbl.index.duplicated()] = [i+'.1' for i in symEnsmbl.loc[symEnsmbl.index.duplicated()].index]
symEnsmbl.index = pd.Index(tmp1)

for set1 in datasets:
    common = [i for i in datasets[set1].var_names if i in symEnsmbl[0]]
    datasets[set1]= datasets[set1][:,common]
    datasets[set1].var_names = pd.Index(symEnsmbl.loc[datasets[set1].var_names,0], name='Ensembl')
    datasets[set1].var_names_make_unique()

# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv(resdir2+'/Whitfield/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
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

# Subset to marker genes
mg1 = []
for set1 in datasets:
    if set1 == 'U5':
        tmp1 = pd.read_csv('markerGenes/'+set1+'_scTransform_Markers_together_030423.csv', header = 0, index_col = 0)
    tmp1 = tmp1.loc[tmp1['avg_log2FC'].ge(0.25),:]
    tmp1 = tmp1[tmp1['cluster'] != 'G1/other']
    tmp1 = tmp1.loc[tmp1['p_val_adj'].le(0.05),:]
    full_filtered_genes = tmp1
    tmp3 = set(list(tmp1['gene']))
    tmp2 = []
    missed = []
    for j in tmp3:
        if j in hgncEnsmbl.index:
            tmp2.append(hgncEnsmbl.loc[j,'Ensembl ID(supplied by Ensembl)'])
        elif j in hgncPrevEnsmbl:
            tmp2.append(hgncPrevEnsmbl[j])
        elif j in hgncAliasEnsmbl:
            tmp2.append(hgncAliasEnsmbl[j])
        else:
            missed.append(j)
    print(set1+' marker genes: '+str(len(tmp2)))
    mg1 += tmp2
    mg1 = list(set(mg1).intersection(datasets[set1].var_names))
    print('U5 marker genes: '+str(len(mg1)))
    # Subset data to marker genes
    datasets[set1] = datasets[set1][:,mg1]


#mgenes = pd.read_csv('ccAFv2_py/ccAF/ccAFv2_genes.csv', index_col = 0)
#merged = pd.merge(full_filtered_genes, mgenes, left_on = 'gene', right_on = 'human_symbol')
#merged.to_csv('ccNN/mgenes_by_cluster.csv')

#------------------------
# Cross validation
#-----------------------

# Parameters for CV
numSamples = 300
nfolds = 10
ncores = 10

# Initialize helper vars/indices for subsetting data (train/test)
nCells = {}
allInds = {}
for set1 in datasets:
    nCells[set1] = datasets[set1].shape[0]
    allInds[set1] = np.arange(0, nCells[set1])

os.makedirs(resdir5+'/classification_method_comparison_redo_with_10_perc_holdout')
resdir7 = resdir5+'/classification_method_comparison_redo_with_10_perc_holdout'

# Precompute testing and training data sets
trainInds = {}
truelab = {}
testInds = {}
mapMe_SVMrejRF = {}
mapMe_KNN = {}
mapMe_ACTINN = {}
mapMe_ccNN = {}
errorACTINN = {}
for set1 in datasets:
    trainInds[set1] = []
    truelab[set1] = []
    testInds[set1] = []
    mapMe_SVMrejRF[set1] = []
    mapMe_KNN[set1] = []
    mapMe_ACTINN[set1] = []
    errorACTINN[set1] = []
    for k in range(nfolds):
        samp1 = sample_without_replacement(nCells[set1], numSamples, random_state = 1234 + k)
        testInds[set1].append(samp1)
        truelab[set1].extend(datasets_mg1[set1][testInds[set1][k],:].obs['new_clusters'])
        trainInds[set1].append(np.setdiff1d(allInds[set1], samp1))
        mapMe_SVMrejRF[set1].append([k, datasets_mg1[set1][trainInds[set1][k],:], datasets_mg1[set1].obs['new_clusters'][trainInds[set1][k]], datasets_mg1[set1][testInds[set1][k],:]])
        mapMe_KNN[set1].append([k, np.setdiff1d(allInds, samp1), samp1, datasets_mg1[set1]])
        mapMe_ACTINN[set1].append([k, np.setdiff1d(allInds, samp1), samp1, datasets_mg1[set1]])

#################
### SVMrej CV ###
#################

if not exists(resdir7+'/CV_classification_report.csv'):
    # SVMrej multiprocessable function
    def runSVMrej(params):
        print('SVMrej round '+str(params[0])+'...')
        svmRej = cl.Classifier_SVMrej(params[1], params[2])
        testPredLbls = svmRej.predict_labels(params[3])
        return [params[0], testPredLbls]
    # Cross validation within dataset
    print('\nSVMrej cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(runSVMrej, mapMe_SVMrejRF[set1])
    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(resdir7+'/SVMrej_ccAF_CV_results_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets_mg1[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'SVMrej'
    performDF.to_csv(resdir7+'/SVMrej_CV_classification_report_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(resdir7+'/SVMrej_CV_classification_report_mean_'+str(datasets_mg1[set1]._n_vars)+'_'+str(k)+'_'+set1+'.csv')

#############
### RF CV ###
#############

if not exists(resdir7+'/CV_classification_report.csv'):
    # RF multiprocessable function
    def runRF(params):
        print('RF round '+str(params[0])+'...')
        RF = cl.Classifier_RF(params[1], params[2])
        testPredLbls = RF.predict_labels(params[3])
        return [params[0], testPredLbls]
    # Cross validation within dataset
    print('\nRF cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(runRF, mapMe_SVMrejRF[set1])
    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(resdir7+'/RFpy_ccAF_CV_results_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets_mg1[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'RFpy'
    performDF.to_csv(resdir7+'/RFpy_CV_classification_report_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(resdir7+'/RFpy_CV_classification_report_mean_'+str(datasets_mg1[set1]._n_vars)+'_'+str(k)+'_'+set1+'.csv')


##############
### KNN CV ###
##############

if not exists(resdir7+'/CV_classification_report.csv'):
    # Cross validation within dataset
    pred = []
    print('\nKNN cross-validation...')
    for k in range(nfolds):
        print('KNN round '+str(k)+'...')
        KNN = cl.Classifier_KNN(datasets_mg1[set1][trainInds[set1][k],:], 'new_clusters')
        testPredLbls = KNN.predict_labels(datasets_mg1[set1][testInds[set1][k],:])
        pred.extend(testPredLbls)
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(resdir7+'/KNN_ccAF_CV_results_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets_mg1[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'KNN'
    performDF.to_csv(resdir7+'/KNN_CV_classification_report_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    #errorccNN.append(DF['Equal'].value_counts(normalize=True))
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(resdir7+'/KNN_CV_classification_report_mean_'+str(datasets_mg1[set1]._n_vars)+'_'+str(k)+'_'+set1+'.csv')


#######################################################################
### ACTINN CV - need to run on different docker / use classifieresV3 ##
#######################################################################

if not exists(resdir7+'/CV_classification_report.csv'):
    # Cross validation within dataset
    pred = []
    print('\nACTINN cross-validation...')
    for k in range(nfolds):
        print('ACTINN round '+str(k)+'...')
        ACTINN = cl.Classifier_ACTINN(train = datasets_mg1[set1][trainInds[set1][k],:], label = 'new_clusters')
        testPredLbls = ACTINN.predict_labels(datasets_mg1[set1][testInds[set1][k],:])
        pred.extend(testPredLbls)
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(resdir7+'/ACTINN_ccAF_CV_results_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets_mg1[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'ACTINN'
    performDF.to_csv(resdir7+'/ACTINN_CV_classification_report_'+str(datasets_mg1[set1]._n_vars)+'.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(resdir7+'/ACTINN_CV_classification_report_mean_'+str(datasets_mg1[set1]._n_vars)+'_'+str(k)+'_'+set1+'.csv')
