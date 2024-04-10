##########################################################
## ccAFv2:  CV compare ccAF on U5s                      ##
##          with ccSeurat (Phase) as reference          ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/:/files' cplaisier/scrna_seq_velocity

import numpy as np
import pandas as pd
import scanpy as sc
import ccAF
from sklearn.utils.random import sample_without_replacement
import os

def flatten(xss):
    return [x for xs in xss for x in xs]

# Set parameters for analysis
nfolds = 10
ncores = 10

# Load data
data1 = sc.read_loom('data/filtered/final/U5_filtered_ensembl.loom')
data1.obs['Cell Labels'] = data1.obs['ccAF']
seurat_calls = pd.read_csv('U5_ccseurat_calls_101023.csv', index_col = 0)
data1.obs['Phase'] = seurat_calls['x']

# Initialize helper vars/indices for subsetting data (test)
barcodes = data1.obs_names
nCells = len(data1)
numSamples = round(0.8*nCells)
allInds = np.arange(0, nCells)
truelab = []
predlab = []
for k in range(nfolds):
    samp1 = sample_without_replacement(nCells, numSamples, random_state = 1234 + k)
    data_sub = data1[samp1,:]
    truelab.append(list(data_sub.obs['Phase']))
    data_sub.obs['ccAF'] = ccAF.ccAF.predict_labels(data_sub)
    predlab.append(list(data_sub.obs['ccAF']))


os.makedirs('cross_validation/ccaf')
DF = pd.DataFrame({'True Labels':flatten(truelab), 'ccAF_pred':flatten(predlab)})
DF.to_csv('cross_validation/ccaf/CV_classification_report_020824_with_Phase_as_ref.csv')

# Just reclassify with all cells
data1 = sc.read_loom('data/filtered/final/U5_filtered_ensembl.loom')
data1.obs['Cell Labels'] = data1.obs['ccAF']
data1.obs['ccAF_pred'] = ccAF.ccAF.predict_labels(data1)
data1.obs['ccAF_pred'].to_csv('U5_ccAF_calls_020224.csv')
