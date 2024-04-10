##########################################################
## ccAFv2:  CV compare cell cycle classifiers on BT322  ##
##          with Phase as ref metrics and plotting      ##
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

# docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccnn

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap

#  metrics
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
from scipy import stats

##################
# Set up section #
##################

## Load datasets and concatenate
tags = ['ccafv2', 'tricycle', 'recat', 'ccschwabe', 'peco', 'cyclone']
resdir = 'cross_validation'
datasets = {}
for tag1 in tags:
    datasets[tag1] = pd.read_csv(resdir+'/BT322/'+tag1+'/CV_classification_report_020824_with_Phase_as_ref.csv', index_col = 0)
    datasets[tag1]['classifier'] = tag1
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[0]].name: 'truelab'}, inplace=True)
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[1]].name: 'predlab'}, inplace=True)

df_all = pd.concat([datasets[tags[0]], datasets[tags[1]], datasets[tags[2]], datasets[tags[3]], datasets[tags[4]], datasets[tags[5]]], axis = 0)
# change Unknown to NaN
df_all['predlab'].replace({'Unknown': np.nan, 'G2M': 'G2/M','G2.M':'G2/M', 'M.G1': 'M/Early G1'}, inplace=True)
#df_all.dropna(inplace = True)


# Tabulate AMIs
ami_scores = {}
for tag1 in tags:
    defined_cell_states = df_all[df_all['classifier'] == tag1]['truelab']
    compare_cell_states = df_all[df_all['classifier'] == tag1]['predlab']
    ami_scores[tag1] = adjusted_mutual_info_score(defined_cell_states[compare_cell_states.dropna().index], compare_cell_states.dropna())


lst = [7, 4, 6, 5, 4, 3, 8]
x = np.array(lst)
y = np.array(list(ami_scores.values()))
y = list(y)
y.append(ccaf_ami)

m = (y[1]-y[2])/(x[1]-x[2])
b = y[1]-(m*x[1])

descrip = np.array(['ccAFv2', 'tricycle', 'recat', 'ccSchwabe', 'peco', 'cyclone', 'ccAF'])
df2 = pd.DataFrame(y, descrip)
df2.rename(columns={0:'AMI'}, inplace=True)
df2.to_csv(resdir+'/BT322_AMI_seurat_reference.csv')

# Plot
fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(6,6))
ax.errorbar(x, y,
            #xerr=xerr1,
            #yerr=yerr1,
            fmt='o', ecolor = 'black')
ax.set_xlabel('# of classes')
ax.set_ylabel('Adjusted Mutual Information Score')
ax.axline((0, b), slope=m, color='C0', label='by slope')
ax.set_ylim(0,1)
ax.set_xlim(0,8)
for i, txt in enumerate(descrip):
    ax.annotate(txt, (x[i], y[i]))
#matplotlib.pyplot.errorbar(np.array(x),np.array(y), np.array(xerr), np.array(yerr), 'o')
#plt.errorbar(x, y, yerr=yerr1, fmt='o')
plt.savefig(resdir+'/BT322_AMI_no_classes_together_seurat_ref_021424.pdf')
plt.clf()







# long way
"""
# Run comparison tests
results = {}
cells_predicted = {}
defined_cell_states = df_all['truelab']
compare_cell_states = df_all['predlab']
numSamples = round(2962*0.8) # removed G1/other cells
nfolds = 10
for tag1 in tags:
    results[tag1] = []
    cells_predicted[tag1] = []
    defined_cell_states = df_all[df_all['classifier'] == tag1]['truelab']
    compare_cell_states = df_all[df_all['classifier'] == tag1]['predlab']
    for k in range(nfolds):
        bind = int(numSamples*k)
        eind = int(numSamples*(k+1))
        truelab = defined_cell_states.iloc[bind:eind]
        #truelab = compare_cell_states.iloc[bind:eind]
        predlab = compare_cell_states.iloc[bind:eind]
        # Adjusted mutual score; drop unknowns
        results[tag1].append(adjusted_mutual_info_score(truelab[predlab.dropna().index], predlab.dropna()))
        # Cells predicted
        cells_predicted[tag1].append(1 - sum(predlab.isna())/len(predlab))

col1 = []
col2 = []
col3 = []
for tag1 in tags:
    tmp1 = pd.DataFrame(range(nfolds), results[tag1]).reset_index().rename(columns = {'index':'adjusted_mutual_score'})
    tmp2 = pd.DataFrame(range(nfolds), cells_predicted[tag1]).reset_index().rename(columns = {'index':'cells_predicted'})
    df = pd.concat([tmp1['adjusted_mutual_score'], tmp2['cells_predicted']], axis=1)
    df['classifier'] = tag1
    df.reset_index(inplace=True)
    col1.append(list(df['adjusted_mutual_score']))
    col2.append(list(df['cells_predicted']))
    col3.append(list(df['classifier']))

col1 = [item for sublist in col1 for item in sublist]
col2 = [item for sublist in col2 for item in sublist]
col3 = [item for sublist in col3 for item in sublist]

df2 = pd.DataFrame([col1, col2, col3]).T
df2.rename(columns={0:'adjusted_mutual_score', 1:'cells_predicted', 2: 'classifier'}, inplace=True)
#df2['data'] = data1
for metric1 in ['adjusted_mutual_score', 'cells_predicted']:
    sns.set(style="whitegrid", font_scale=3)
    fig, ax1 = plt.subplots(figsize=(40,20))
    sns.boxplot(hue = 'classifier', y = metric1, x = "data", data = df2, ax = ax1, palette = "husl")
    ax1.set(ylabel=metric1)
    ax1.set_ylim(0,1)
    #plt.savefig(resdir+'/'+data1+'_'+metric1+'_classifier_comparison_111623.pdf')
    #plt.savefig(resdir+'/'+data1+'_'+metric1+'_classifier_comparison_010424.pdf')
    plt.savefig(resdir+'/'+data1+'_'+metric1+'_classifier_comparison_020524.pdf')
    plt.clf()

ccaf_pred = pd.read_csv(resdir+'/BT322/ccaf/CV_classification_report_020824_with_Phase_as_ref.csv', index_col = 0)
ccaf_pred['classifier'] = tag1
ccaf_pred.rename(columns = {ccaf_pred[ccaf_pred.columns[0]].name: 'truelab'}, inplace=True)
ccaf_pred.rename(columns = {ccaf_pred[ccaf_pred.columns[1]].name: 'predlab'}, inplace=True)
truelab = ccaf_pred['truelab']
predlab = ccaf_pred['predlab']
ccaf_ami = adjusted_mutual_info_score(truelab, predlab)


median_ami = []
median_cp = []
for classifier1 in list(df2['classifier'].unique()):
    median_ami.append(np.median(df2[df2['classifier'] == classifier1]['adjusted_mutual_score']))
    median_cp.append(np.median(df2[df2['classifier'] == classifier1]['cells_predicted']))

median_cp_perc = [i * 100 for i in median_cp]

x = np.array(median_cp_perc)
y = np.array(median_ami)
#xerr1 = np.array([0.001, 0, 0.002, 0, 0.005, 0.003, 0])
#yerr1 = np.array([0.005, 0.005, 0.006, 0.006, 0.004, 0.004, 0.004])
#descrip = np.array(['ccAF', 'ccAFv2', 'ccSeurat', 'tricycle', 'recat', 'ccSchwabe', 'peco', 'cyclone'])
descrip = np.array(['ccAFv2', 'tricycle', 'recat', 'ccSchwabe', 'peco', 'cyclone'])
fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10,10))
ax.errorbar(x, y,
            #xerr=xerr1,
            #yerr=yerr1,
            fmt='o', ecolor = 'black')
ax.set_xlabel('Cells Predicted (%)')
ax.set_ylabel('Adjusted Mutual Information Score')
ax.set_ylim(0,1.1)
ax.set_xlim(50,100)
for i, txt in enumerate(descrip):
    ax.annotate(txt, (x[i], y[i]))
#matplotlib.pyplot.errorbar(np.array(x),np.array(y), np.array(xerr), np.array(yerr), 'o')
#plt.errorbar(x, y, yerr=yerr1, fmt='o')
plt.savefig(resdir+'/BT322/AMI_cellsPredicted_together_020824.pdf')
plt.clf()
"""
