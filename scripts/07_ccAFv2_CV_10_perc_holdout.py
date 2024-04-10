##########################################################
## ccAFv2: ccAFv2 CV 10 percent holdout                 ##
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

###############
## Imports   ##
###############
# General
import numpy as np
import pandas as pd
import scanpy as sc
import os
from scipy.sparse import isspmatrix
import pickle
import calendar
import time
time_stamp = calendar.timegm(time.gmtime())

# ccAFv2
from sklearn.preprocessing import StandardScaler, LabelEncoder
import numpy as np
import tensorflow as tf
from tensorflow import keras

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
from sklearn.metrics import adjusted_mutual_info_score

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

#################
## Classifiers ##
#################

class Classifier_ccAFv2:
    """A class designed to facilitate our new ccNN classifier construction and use."""
    def __init__(self, training_sets, label_sets, variable_genes = 1000, epochs = 10, training_iterations = 10, validation_split = 0.2, activation = 'relu', dropout_rate = 0.4, model_specs = [200,100], model_name = calendar.timegm(time.gmtime())):
        self.epochs = epochs
        self.training_iterations = training_iterations
        self.validation_split = validation_split
        self.activation = activation
        self.dropout_rate = dropout_rate
        self.model_name = model_name
        self.classifier, self.label_encoder, self.genes = self.__build_classifier(training_sets, label_sets, variable_genes, epochs, training_iterations, validation_split, activation, dropout_rate, model_specs)
    # Prepare data
    def __prep_data_sets(self, data_sets):
        data_sets_scaled = {}
        scaler = StandardScaler()
        for data in data_sets:
            # Make indicies unique
            data_sets[data].var_names_make_unique()
            # Remove all genes with zero counts
            sc.pp.filter_genes(data_sets[data], min_cells=1)
            if isspmatrix(data_sets[data].X):
                data_sets_scaled[data] = pd.DataFrame(data_sets[data].X.todense(), index = data_sets[data].obs_names, columns = data_sets[data].var_names)
            else:
                data_sets_scaled[data] = pd.DataFrame(data_sets[data].X, index = data_sets[data].obs_names, columns = data_sets[data].var_names)
            data_sets_scaled[data] = _scale(data_sets_scaled[data])
        return data_sets_scaled
    # Prepare test data for predicting
    def __prep_predict_data(self, data):
        # Remove all genes with zero counts
        data.var_names_make_unique()
        sc.pp.filter_genes(data, min_cells=1)
        # Restrict to classifier genes
        data2 = data[:,list(set(self.genes).intersection(data.var_names))]
        # Scale data
        scaler = StandardScaler()
        if isspmatrix(data.X):
            data2 = pd.DataFrame(data2.X.todense(), index = data2.obs_names, columns = data2.var_names)
        else:
            data2 = pd.DataFrame(data2.X, index = data2.obs_names, columns = data2.var_names)
        data3 = pd.DataFrame(_scale(data2), index = data2.index, columns = data2.columns)
        # Add minimum values for missing genes
        missing = set(self.genes).difference(data3.columns)
        if len(missing)>0:
            data4 = pd.concat([data3, pd.DataFrame(data3.values.min(), index=data3.index, columns = missing)], axis=1)
            return data4[list(self.genes)]
        else:
            return data3
    # Build classifier
    def __build_classifier(self, training_sets, label_sets, variable_genes, epochs, training_iterations, validation_split, activation, dropout_rate, model_specs):
        # Get common genes
        gene_sets = [set(training_sets[set1].var_names) for set1 in training_sets]
        common_genes = list(gene_sets[0].intersection(*gene_sets))
        training_sets_common = {set1:training_sets[set1][:,common_genes] for set1 in training_sets}
        # Subset to highly variable genes
        training_sets_hvgs = training_sets_common
        # Make label encoder and encode labels
        label_encoder = LabelEncoder()
        label_encoder.fit([j for label_set in label_sets for j in label_sets[label_set]])
        label_sets_prepared = {set1:label_encoder.transform(label_sets[set1]) for set1 in label_sets}
        # Convert into pandas DataFrame
        training_sets_prepared = self.__prep_data_sets(training_sets_hvgs)
        classifier = self.__train_model(training_sets_prepared, label_sets_prepared, training_sets_hvgs[list(training_sets_hvgs.keys())[0]].shape[1], epochs, training_iterations, validation_split, activation, dropout_rate, model_specs)
        return classifier, label_encoder, common_genes
    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff=0.5):
        pred_data = self.__prep_predict_data(new_data)
        probabilities = self.__predict_new_data(pred_data, self.classifier)
        labels = self.label_encoder.inverse_transform([np.argmax(i) for i in probabilities])
        labels[np.where([np.max(i) < cutoff for i in probabilities])] = np.nan
        return labels, probabilities
    # Save our weights from model
    def save_weights(self, file_name):
        self.classifier.save_weights(file_name)
    # Save entire model as an hdf5 file
    def save_model_hdf5(self, file_name):
        self.classifier.save_model_hdf5(file_name)
    # Train the ccAFv2 model
    def __train_model(self, training_sets, label_sets, input_width, epochs = 10, training_iterations = 5, validation_split = 0.2, activation = 'relu', dropout_rate = 0.4, model_specs = [200,100]):
        # Build Neural Network model
        model = keras.models.Sequential()
        model.add(keras.layers.Dense(input_width, input_dim=input_width, activation=activation))
        for layer1 in range(len(model_specs)):
            model.add(keras.layers.Dense(model_specs[layer1], activation=activation))
            model.add(keras.layers.Dropout(dropout_rate))
        model.add(keras.layers.Dense(len(set([j for i in label_sets for j in label_sets[i]])), activation='softmax'))
        # Compile model -> sgd is best optimizer that was tested
        model.compile(loss=tf.keras.losses.CategoricalCrossentropy(), optimizer='sgd', metrics=['accuracy'])
        # Train NN five times through per dataset using online learning with 10 epochs
        for i in range(training_iterations):
            for set1 in training_sets:
                print('Training round '+str(i+1)+'; Training set: '+str(set1))
                model.fit(np.array(training_sets[set1]), tf.keras.utils.to_categorical(label_sets[set1]), epochs=epochs, validation_split=validation_split)
        return model
    # Predict ccAFv2 labels for new data
    def __predict_new_data(self, new_data, model):
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
        return model.predict(new_data)

##########################################
## ccAFv2 CV 10% hold out
##########################################
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

# Parameters
numSamples = 300
nfolds = 10
ncores = 10
layers = ['600_200']
if layers:
    for lay1 in layers:
        print(lay1)
        model_name = 'layers_'+lay1
        model_specs1 = [eval(i) for i in model_name.split('layers_')[1].split('_')]
        print('\nBuilding ccAFv2...')
        ccAFv2 = Classifier_ccAFv2(training_sets = {set1:datasets[set1]}, label_sets = {set1:datasets[set1].obs['new_clusters']}, dropout_rate=0.5, model_specs = model_specs1)
        #------------------------
        # Cross validation
        #-----------------------
        # Initialize for cross validation
        nCells = {}
        allInds = {}
        trainInds = {}
        testInds = {}
        truelab = {}
        pred = {}
        likelihoods = {}
        training_sets = {}
        label_sets = {}
        for set1 in datasets:
            nCells[set1] = datasets[set1].shape[0]
            allInds[set1] = np.arange(0, nCells[set1])
            trainInds[set1] = []
            testInds[set1] = []
            truelab[set1] = []
            pred[set1] = []
            likelihoods[set1] = []
            for k in range(nfolds):
                samp1 = sample_without_replacement(nCells[set1], numSamples, random_state = 1234 + k)
                testInds[set1].append(samp1)
                truelab[set1].extend(datasets[set1][testInds[set1][k],:].obs['new_clusters'])
                trainInds[set1].append(np.setdiff1d(allInds[set1], samp1))
                training_sets[set1] = datasets[set1][trainInds[set1][k],:]
                label_sets[set1] = datasets[set1][trainInds[set1][k],:].obs['new_clusters']
        #################
        ### ccAFv2 CV ###
        #################
        ## Load up data
        if not os.path.exists(resdir2+'/ccAFv2'):
            os.makedirs(resdir2+'/ccAFv2')
        savedir = resdir2+'/ccAFv2'
        if not exists(savedir+'/CV_classification_report.csv'):
            print('\nccAFv2 cross-validation...')
            for k in range(nfolds):
                print('ccAFv2 round '+str(k)+'...')
                for set1 in sets:
                    testPredLbls = ccAFv2.predict_labels(datasets[set1][testInds[set1][k],:])
                    d1 = ccAFv2.predict_labels(datasets[set1])
                    df1 = pd.concat([pd.Series(d1[0]),pd.DataFrame(d1[1])],axis=1)
                    df1.index = datasets[set1].obs_names
                    df1.columns = ['Prediction']+list(ccNN.label_encoder.classes_)
                    pred[set1].extend(testPredLbls[0])
                    likelihoods[set1].extend(testPredLbls[1])
            # Dataframe of true labels, predictions, probabilities for all iterations
            for set1 in sets:
                DF = pd.DataFrame({'True Labels':truelab[set1], 'Predictions':pred[set1]})
                DF = DF[DF['Predictions'] != 'nan']
                DF.to_csv(savedir+'/ccAFv2_CV_results_'+str(datasets[set1]._n_vars)+'_'+set1+'_'+model_name+'_10_perc_holdout.csv')
                # Get classification report for each iteration
                performanceResults = []
                for k in range(nfolds):
                    slice1 = slice((numSamples)*k, (numSamples)*(k+1), 1)
                    true1 = truelab[set1][slice1]
                    pred1 = pred[set1][slice1]
                    s1 = [i for i in range(len(pred1)) if pred1[i] != 'nan']
                    true2 = [true1[i] for i in s1]
                    pred2 = [pred1[i] for i in s1]
                    performanceResults.append(classification_report(true2, pred2, output_dict=True, zero_division=0))
                # Convert into a dataframe
                performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
                states1 = list(pd.DataFrame(datasets[set1].obs['new_clusters'].value_counts()).index)
                performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
                performDF['Classifier'] = 'ccAFv2'
                performDF.to_csv(savedir+'/ccAFv2_CV_classification_report_'+str(datasets[set1]._n_vars)+'_'+set1+'_'+model_name+'_10_perc_holdout.csv')
                comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
                DF["Equal"] = comparison_column
                performDF.index.name = 'index1'
                performDF.groupby(by='index1').mean().to_csv(savedir+'/ccAFv2_CV_classification_report_mean_'+str(datasets[set1]._n_vars)+'_'+set1+'_'+model_name+'_10_perc_holdout.csv')
