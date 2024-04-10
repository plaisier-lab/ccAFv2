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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccnn
#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/scrna_seq_velocity # ACTINN
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
import classifiersV3 as cl
#import classifiersV6clp as cl

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

# Set up
sets = ['U5']
resdir = 'data/'+sets[0]
resdir2 = 'results'
p_cutoffs = {'U5':0.05}

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
        truelab[set1].extend(datasets[set1][testInds[set1][k],:].obs['new_clusters'])
        trainInds[set1].append(np.setdiff1d(allInds[set1], samp1))
        mapMe_SVMrejRF[set1].append([k, datasets[set1][trainInds[set1][k],:], datasets[set1].obs['new_clusters'][trainInds[set1][k]], datasets[set1][testInds[set1][k],:]])
        mapMe_KNN[set1].append([k, np.setdiff1d(allInds, samp1), samp1, datasets[set1]])
        mapMe_ACTINN[set1].append([k, np.setdiff1d(allInds, samp1), samp1, datasets[set1]])

#################
### SVMrej CV ###
#################

if not os.path.exists(resdir2+'/SVMrej'):
    os.makedirs(resdir2+'/SVMrej')

savedir = resdir2+'/SVMrej'
if not exists(savedir+'/CV_classification_report.csv'):
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
    DF.to_csv(savedir+'/SVMrej_CV_results_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'SVMrej'
    performDF.to_csv(savedir+'/SVMrej_CV_classification_report_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(savedir+'/SVMrej_CV_classification_report_mean_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')

#############
### RF CV ###
#############

if not os.path.exists(resdir2+'/RFpy'):
    os.makedirs(resdir2+'/RFpy')

savedir = resdir2+'/RFpy'
if not exists(savedir+'/CV_classification_report.csv'):
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
    DF.to_csv(savedir+'/RFpy_ccAF_CV_results_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'RFpy'
    performDF.to_csv(savedir+'/RFpy_CV_classification_report_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(savedir+'/RFpy_CV_classification_report_mean_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')


##############
### KNN CV ###
##############
if not os.path.exists(resdir2+'/KNN'):
    os.makedirs(resdir2+'/KNN')

savedir = resdir2+'/KNN'
if not exists(savedir+'/CV_classification_report.csv'):
    # Cross validation within dataset
    pred = []
    print('\nKNN cross-validation...')
    for k in range(nfolds):
        print('KNN round '+str(k)+'...')
        KNN = cl.Classifier_KNN(datasets[set1][trainInds[set1][k],:], 'new_clusters')
        testPredLbls = KNN.predict_labels(datasets[set1][testInds[set1][k],:])
        pred.extend(testPredLbls)
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(savedir+'/KNN_ccAF_CV_results_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'KNN'
    performDF.to_csv(savedir+'/KNN_CV_classification_report_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    #errorccNN.append(DF['Equal'].value_counts(normalize=True))
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(savedir+'/KNN_CV_classification_report_mean_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')


####################################################################################
### ACTINN CV - need to run on different docker / use classifier classifiersV6clp ##
###################################################################################
"""
if not os.path.exists(resdir2+'/ACTINN'):
    os.makedirs(resdir2+'/ACTINN')

savedir = resdir2+'/ACTINN'
if not exists(savedir+'/CV_classification_report.csv'):
    # Cross validation within dataset
    pred = []
    print('\nACTINN cross-validation...')
    for k in range(nfolds):
        print('ACTINN round '+str(k)+'...')
        ACTINN = cl.Classifier_ACTINN(train = datasets[set1][trainInds[set1][k],:], label = 'new_clusters')
        testPredLbls = ACTINN.predict_labels(datasets[set1][testInds[set1][k],:])
        pred.extend(testPredLbls)
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred})
    DF.to_csv(savedir+'/ACTINN_ccAF_CV_results_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[set1][slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = list(pd.DataFrame(datasets[set1].obs['new_clusters'].value_counts()).index)
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'ACTINN'
    performDF.to_csv(savedir+'/ACTINN_CV_classification_report_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    performDF.index.name = 'index1'
    performDF.groupby(by='index1').mean().to_csv(savedir+'/ACTINN_CV_classification_report_mean_'+str(datasets[set1]._n_vars)+'_'+set1+'_10_perc_holdout.csv')
"""
