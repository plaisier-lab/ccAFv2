##########################################################
## ccAFv2:  CV compare cell cycle classifiers on U5s    ##
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
data1 = 'U5'
tags = ['ccaf', 'ccafv2', 'tricycle', 'recat', 'ccschwabe', 'peco', 'cyclone']
resdir = 'cross_validation'
datasets = {}
for tag1 in tags:
    datasets[tag1] = pd.read_csv(resdir+'/'+tag1+'/CV_classification_report_020824_with_Phase_as_ref.csv')
    datasets[tag1]['classifier'] = tag1
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[1]].name: 'truelab'}, inplace=True)
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[2]].name: 'predlab'}, inplace=True)

df_all = pd.concat([datasets[tags[0]], datasets[tags[1]], datasets[tags[2]], datasets[tags[3]], datasets[tags[4]], datasets[tags[5]], datasets[tags[6]]], axis = 0)
# change Unknown to NaN
df_all['predlab'].replace({'Unknown': np.nan, 'G2M': 'G2/M','G2.M':'G2/M', 'M.G1': 'M/Early G1'}, inplace=True)
#df_all.dropna(inplace = True)

# QUICK WAY - medians
ami_scores = {}
for tag1 in tags:
    defined_cell_states = df_all[df_all['classifier'] == tag1]['truelab']
    compare_cell_states = df_all[df_all['classifier'] == tag1]['predlab']
    ami_scores[tag1] = adjusted_mutual_info_score(defined_cell_states[compare_cell_states.dropna().index], compare_cell_states.dropna())

df = pd.DataFrame.from_dict(ami_scores, orient = 'index')
lst = [8, 7, 4, 6, 5, 4, 3, 7]
x = np.array(lst)
y = np.array(list(ami_scores.values()))
y = list(y)
y.append(0.34168950812765647)

m = (y[2]-y[3])/(x[2]-x[3])
b = y[2]-(m*x[2])

descrip = np.array(['ccAF', 'ccAFv2', 'tricycle', 'recat', 'ccSchwabe', 'peco', 'cyclone', 'cell labels'])
df2 = pd.DataFrame(y, descrip)
df2.rename(columns={0:'AMI'}, inplace=True)
df2.to_csv(resdir+'/U5_AMI_seurat_reference.csv')


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
plt.savefig(resdir+'/U5_AMI_no_classes_together_seurat_ref_021424.pdf')
plt.clf()





"""
# Try different metrics for comparison
# Load in U5 SCT expression data
data = pd.read_csv('U5_SCT_expression_861_genes.csv', index_col=0)
#data = pd.read_csv('U5_SCT_data_expression_861_genes.csv', index_col=0)

datasets = {}
tags = ['cell_labels','ccAFv2', 'ccAF', 'ccseurat', 'tricycle', 'recat', 'schwabe', 'peco', 'cyclone']
for tag1 in tags:
    if tag1 in ['cell_labels', 'tricycle', 'peco', 'ccAFv2', 'ccAF']:
        datasets[tag1] = pd.read_csv('U5_'+tag1+'_calls_020224.csv')
    else:
        datasets[tag1] = pd.read_csv('U5_'+tag1+'_calls_101023.csv')
    datasets[tag1]['classifier'] = tag1
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[1]].name: 'predlab'}, inplace=True)

df_all = pd.concat([datasets[tags[0]], datasets[tags[1]], datasets[tags[2]], datasets[tags[3]], datasets[tags[4]], datasets[tags[5]], datasets[tags[6]], datasets[tags[7]], datasets[tags[8]]], axis = 0)
# change Unknown to NaN
df_all['predlab'].replace({'Unknown': np.nan, 'G2M': 'G2/M','G2.M':'G2/M', 'M.G1': 'M/Early G1'}, inplace=True)
#df_all['predlab'].replace({'G2M': 'G2/M','G2.M':'G2/M', 'M.G1': 'M/Early G1'}, inplace=True)


data2 = data.T
data2.reset_index(inplace = True)
scores = {}
for tag1 in tags:
    tmp = df_all[df_all['classifier'] == tag1]['predlab']
    compare_cell_states = tmp.dropna()
    comp_exprs = data2.loc[compare_cell_states.dropna().index]
    comp_exprs.set_index('index', inplace = True)
    #scores[tag1] = sklearn.metrics.davies_bouldin_score(comp_exprs, compare_cell_states)
    #scores[tag1] = sklearn.metrics.davies_bouldin_score(data2, tmp)
    scores[tag1] = sklearn.metrics.calinski_harabasz_score(comp_exprs, compare_cell_states)
    #scores[tag1] = sklearn.metrics.calinski_harabasz_score(data2, tmp)


df = pd.DataFrame.from_dict(scores, orient = 'index')
#df2 = df.rename(columns={0:'davies_bouldin_score'})
#df2 = df.rename(columns={0:'calinski_harabasz_score'})

lst = [7, 7, 8, 3, 4, 6, 5, 4, 3]
x = np.array(lst)
y = np.array(list(scores.values()))
descrip = np.array(tags)
fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10,10))
ax.errorbar(x, y,
            #xerr=xerr1,
            #yerr=yerr1,
            fmt='o', ecolor = 'black')
ax.set_xlabel('# of classes')
ax.set_ylabel('calinski_harabasz_score')
ax.set_ylim(0,200)
ax.set_xlim(0,8)
for i, txt in enumerate(descrip):
    ax.annotate(txt, (x[i], y[i]))
#matplotlib.pyplot.errorbar(np.array(x),np.array(y), np.array(xerr), np.array(yerr), 'o')
#plt.errorbar(x, y, yerr=yerr1, fmt='o')
plt.savefig('cross_validation/U5_CHS_num_classes_together_021424.pdf')
plt.clf()

"""
