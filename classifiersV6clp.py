##########################################################
## OncoMerge:  classifiers.py                           ##
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

##########################################
## Load Python packages for classifiers ##
##########################################

# General
import numpy as np
import pandas as pd
import os
from scipy.sparse import isspmatrix

# sklearn: Support Vector Machine (SVM) rejection
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import BaggingClassifier

# sklearn: Random Forest (RF)
from sklearn.ensemble import RandomForestClassifier

# ACTINN
#import actinn

# ccNN
from sklearn.preprocessing import StandardScaler, LabelEncoder
import ccNN

# scanpy: K-Nearest Neighbors (KNN)
import scanpy as sc
import anndata


###############
## Functions ##
###############
def scale(df1):
    return (df1.subtract(df1.mean(axis=1),axis=0)).div(df1.std(axis=1),axis=0)

#################
## Classifiers ##
#################

class Classifier_SVMrej:
    """A class designed to facilitate SVMrej classifier construction
    and use."""
    def __init__(self, data, labels, cutoff = 0.7):
        self.data = data
        self.labels = labels
        self.cutoff = cutoff
        self.genes = data.var_names
        self.classifier = self.__build_classifier(data, labels)

    # Build classifier
    def __build_classifier(self, data, labels):
        train = self.__prep_data(data)
        Classifier = LinearSVC(max_iter = 100000)
        clf = CalibratedClassifierCV(Classifier)
        clf.fit(train, labels)
        return clf

    # Predict probability
    def predict_prob(self, newData):
        return np.max(self.classifier.predict_proba(newData), axis=1)

    # Prepare data
    def __prep_data(self, data):
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names)
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if isspmatrix(test_data.X):
            data = pd.DataFrame(test_data.X.todense(), index = test_data.obs_names, columns = test_data.var_names)
        else:
            data = pd.DataFrame(test_data.X, index = test_data.obs_names, columns = test_data.var_names)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(data.values.min(),index=data.index, columns = missing)], axis=1)
        return data[list(self.genes)]

    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff = None):
        pred_data = self.__prep_predict_data(new_data)
        labels = self.classifier.predict(pred_data)
        if cutoff == None and not self.cutoff == None:
            cutoff = self.cutoff
        if not cutoff == None:
            probs = self.predict_prob(pred_data)
            unlabeled = np.where(probs < cutoff)
            labels[unlabeled] = 'Unknown'
        return labels


class Classifier_KNN:
    """A class designed to facilitate scanpy ingetst based
    K-nearest neighbor classifier construction and use."""
    def __init__(self, data, label):
        self.data = data
        self.genes = data.var_names
        self.label = label
        # self.classifier = self.__build_classifier(data, labels)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if len(missing)>0:
            data = pd.concat([pd.DataFrame(test_data.X, index=test_data.obs_names, columns=test_data.var_names), pd.DataFrame(data.values.min(),index=test_data.obs_names, columns = missing)], axis=1)
            data = data[list(self.genes)]
            data_sc = anndata.AnnData(X=data.to_numpy())
            data_sc.var_names = data.columns
            data_sc.obs_names = data.index
            return data_sc
        else:
            return test_data

    # Predict labels
    def predict_labels(self, new_data):
        # Subset based on common gene names
        adata_query = self.__prep_predict_data(new_data)
        var_names = self.data.var_names.intersection(adata_query.var_names)
        adata_ref = self.data[:,var_names]
        adata_query = adata_query[:,var_names]

        # Run embedding anlaysis using subset
        sc.pp.pca(adata_ref)
        sc.pp.neighbors(adata_ref)
        sc.tl.umap(adata_ref)

        # Map the identifiers from the reference dataset to the query dataset
        sc.tl.ingest(adata_query, adata_ref, obs=self.label)

        # Save the results in whitfield data object
        return adata_query.obs[self.label]


class Classifier_RF:
    """A class designed to facilitate RF classifier construction
    and use. Can also be used for """
    def __init__(self, data, labels, cutoff = None):
        self.data = data
        self.genes = data.var_names
        self.labels = labels
        self.cutoff = cutoff
        self.classifier = self.__build_classifier()

    # Build classifier
    def __build_classifier(self):
        train = self.__prep_data(self.data)
        clf = RandomForestClassifier(n_estimators=500, oob_score=True)
        clf.fit(train, self.labels)
        return clf

    # Predict probability
    def predict_prob(self, newData):
        return np.max(self.classifier.predict_proba(newData), axis=1)

    # Prepare data
    def __prep_data(self, data):
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names)
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if isspmatrix(test_data.X):
            data = pd.DataFrame(test_data.X.todense(), index = test_data.obs_names, columns = test_data.var_names)
        else:
            data = pd.DataFrame(test_data.X, index = test_data.obs_names, columns = test_data.var_names)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(data.values.min(),index=data.index, columns = missing)], axis=1)
        return data[list(self.genes)]

    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff = None):
        pred_data = self.__prep_predict_data(new_data)
        labels = self.classifier.predict(pred_data)
        if cutoff == None and not self.cutoff == None:
            cutoff = self.cutoff
        if not cutoff == None:
            probs = self.predict_prob(pred_data)
            unlabeled = np.where(probs < cutoff)
            labels[unlabeled] = 'Unknown'
        return labels



class Classifier_ACTINN:
    """A class designed to facilitate ACTINN classifier construction
    and use. Can also be used for """
    def __init__(self, train, label, learning_rate = 0.0001, num_epochs = 200, minibatch_size = 128, print_cost = True, output_probability = False):
        self.train = train
        self.label = label
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.minibatch_size = minibatch_size
        self.print_cost = print_cost
        self.output_probability = output_probability
        self.label = label
        self.classifier, self.label_to_type_dict, self.genes = self.__build_classifier()

    # Prepare data
    def __prep_data(self, data):
        # Make indicies unique for
        data.var_names_make_unique()
        # Remove all genes with zero counts
        sc.pp.filter_genes(data, min_cells=1)
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names).T
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names).T

    # Prepare test data for predicting
    def __prep_predict_data(self, data):
        missing = set(self.genes).difference(data.index)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(data.values.min(),index=missing, columns = data.columns)])
        return data.loc[list(self.genes)]

    # Build classifier
    def __build_classifier(self):
        train = self.train
        # Convert into pandas DataFrame
        train_data = self.__prep_data(train)
        labels = self.train.obs[self.label]
        clf, label_to_type_dict, genes = actinn.train_model(train_data, labels, learning_rate = self.learning_rate, num_epochs = self.num_epochs, minibatch_size = self.minibatch_size, print_cost = self.print_cost)
        return clf, label_to_type_dict, genes

    # Predict labels with rejection
    def predict_labels(self, newData):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        labels = actinn.predict_new_data(pred_data, self.classifier, self.label_to_type_dict, self.genes)
        return list(labels['celltype'])


class Classifiers_Simple:
    """A class designed to facilitate medioid (argmin) construction
    and use. Can also be used for """
    def __init__(self, data, markergenes, labels, cutoff = None):
        self.data = data
        self.genes = data.var_names
        self.markergenes = markergenes # markergenes is a list
        self.labels = labels
        self.classifier = self.__build_classifier()
    # Prepare data
    def __prep_data(self, data, markergenes):
        # Subset genes to marker genes
        data = data[:,data.var_names[data.var_names.isin(markergenes)]]
        # Remove all genes with zero counts
        sc.pp.filter_cells(data, min_genes=1)
    # Build classifier (argmin medioids)
    def __build_classifier(self):
        # Generate medioids (argmin)
        train = self.train
        train_data = self.__prep_data(train)
        labels = self.train.obs[self.label]
        medioids = pd.DataFrame(index=train_data.var_names, columns=labels.unique())
        for lab1 in list(labels.unique()):
            cells = data[data.obs[self.label]==lab1]
            if len(cells.obs) > 1:
                ind = np.argmin(cells.X.sum(axis=0))
                medioids[lab1] = cells[ind].X.todense().T
        return medioids
    # Prepare test data for predicting
    def __prep_predict_data(self, newData):
        missing = set(self.genes).difference(newData.index)
        if len(missing)>0:
            newData = pd.concat([newData, pd.DataFrame(data.values.min(),index=missing, columns = newData.columns)])
            newData = newData[:,newData.var_names[newData.var_names.isin(medioids.index)]] #subset to marker genes
        return newData.loc[list(self.genes)]
    # Predict labels
    def predict_labels(self, newData):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        corrDF = pd.DataFrame(index=pred_data.obs.index, columns=medioids.columns)
        pred = []
        prob = []
        for cell1 in list(pred_data.obs.index):
            geneExp = pd.DataFrame(index=pred_data.var_names, columns=pred_data[cell1].obs.index)
            geneExp[cell1] = pred_data[cell1].X.todense().T
            tmp = geneExp.reindex(medioids.index)
            for state in medioids.columns:
                corr, _ = spearmanr(tmp,medioids[state])
                corrDF.loc[cell1][state] = corr
            tmp2 = corrDF.loc[cell1]==corrDF.loc[cell1].max()
            pred.append(pd.DataFrame(tmp2[tmp2]).reset_index()['index'][0])
            prob.append(corrDF.loc[cell1][pd.DataFrame(tmp2[tmp2]).reset_index()['index'][0]])
        return pred
        return prob

class Classifier_ccNN:
    """A class designed to facilitate our new ccNN classifier construction and use."""
    def __init__(self, training_sets, label_sets, variable_genes = 1000, epochs = 10, training_iterations = 5, validation_split = 0.2, activation = 'relu', dropout_rate = 0.4):
        self.epochs = epochs
        self.training_iterations = training_iterations
        self.validation_split = validation_split
        self.activation = activation
        self.dropout_rate = dropout_rate
        self.classifier, self.label_encoder, self.genes = self.__build_classifier(training_sets, label_sets, variable_genes, epochs, training_iterations, validation_split, activation, dropout_rate)

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
            #data_sets_scaled[data] = scaler.fit_transform(data_sets_scaled[data])
            data_sets_scaled[data] = scale(data_sets_scaled[data])
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
        #data3 = pd.DataFrame(scaler.fit_transform(data2), index = data2.index, columns = data2.columns)
        data3 = pd.DataFrame(scale(data2), index = data2.index, columns = data2.columns)
        data3.to_csv('training_data_post_standard_scale.csv')
        #data3 = data2
        # Add 0 values for missing genes
        missing = set(self.genes).difference(data3.columns)
        if len(missing)>0:
            data4 = pd.concat([data3, pd.DataFrame(data3.values.min(), index=data3.index, columns = missing)], axis=1)
            return data4[list(self.genes)]
        else:
            return data3

    # Build classifier
    def __build_classifier(self, training_sets, label_sets, variable_genes, epochs, training_iterations, validation_split, activation, dropout_rate):
        # Get common genes
        gene_sets = [set(training_sets[set1].var_names) for set1 in training_sets]
        common_genes = list(gene_sets[0].intersection(*gene_sets))
        training_sets_common = {set1:training_sets[set1][:,common_genes] for set1 in training_sets}
        # Subset to highly variable genes
        #highly_variable_genes = [pd.Series(np.ravel(np.var(training_sets_common[set1].X.todense(),axis=0)), index=training_sets_common[set1].var_names).sort_values(ascending=False).index[range(variable_genes)] for set1 in training_sets_common]
        #hvgs = list(set([gene for set1 in highly_variable_genes for gene in set1]))
        #training_sets_hvgs = {set1:training_sets_common[set1][:,hvgs] for set1 in training_sets_common}
        training_sets_hvgs = training_sets_common
        # Make label encoder and encode labels
        label_encoder = LabelEncoder()
        label_encoder.fit([j for label_set in label_sets for j in label_sets[label_set]])
        label_sets_prepared = {set1:label_encoder.transform(label_sets[set1]) for set1 in label_sets}
        # Convert into pandas DataFrame
        training_sets_prepared = self.__prep_data_sets(training_sets_hvgs)
        classifier = ccNN.train_model(training_sets_prepared, label_sets_prepared, training_sets_hvgs[list(training_sets_hvgs.keys())[0]].shape[1], epochs, training_iterations, validation_split, activation, dropout_rate)
        #classifier = ccNN.train_model(training_sets_prepared, label_sets_prepared, len(hvgs), epochs, training_iterations, validation_split, activation, dropout_rate)
        return classifier, label_encoder, common_genes
        #return classifier, label_encoder, hvgs

    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff=0.5):
        pred_data = self.__prep_predict_data(new_data)
        probabilities = ccNN.predict_new_data(pred_data, self.classifier)
        labels = self.label_encoder.inverse_transform([np.argmax(i) for i in probabilities])
        labels[np.where([np.max(i) < cutoff for i in probabilities])] = np.nan
        return labels, probabilities

    #
    def save_weights(self, file_name):
        self.classifier.save_weights(file_name)

    def save_model_hdf5(self, file_name):
        self.classifier.save_model_hdf5(file_name)
