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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccnn

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
resdir = 'compare_classifiers'
datasets = {}
for tag1 in tags:
    datasets[tag1] = pd.read_csv(resdir+'/'+data1+'/'+tag1+'_CV_classification_report_with_Phase_as_ref.csv', index_col = 0)
    datasets[tag1]['classifier'] = tag1
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[0]].name: 'truelab'}, inplace=True)
    datasets[tag1].rename(columns = {datasets[tag1][datasets[tag1].columns[1]].name: 'predlab'}, inplace=True)

df_all = pd.concat([datasets[tags[0]], datasets[tags[1]], datasets[tags[2]], datasets[tags[3]], datasets[tags[4]], datasets[tags[5]], datasets[tags[6]]], axis = 0)
# change Unknown to NaN
df_all['predlab'].replace({'Unknown': np.nan, 'G2M': 'G2/M','G2.M':'G2/M', 'M.G1': 'M/Early G1'}, inplace=True)
#df_all.dropna(inplace = True)

# Compare adjusted mutual information score for each classifier
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
y.append(0.34168950812765647) # cell labels ami with ref
descrip = np.array(['ccAF', 'ccAFv2', 'tricycle', 'recat', 'ccSchwabe', 'peco', 'cyclone', 'cell labels'])
df2 = pd.DataFrame(y, descrip)
df2.rename(columns={0:'AMI'}, inplace=True)
df2.to_csv(resdir+'/'+data1+'/'+data1+'_AMI_seurat_reference.csv')

fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10,8))
ax.errorbar(x, y,
            fmt='o', ecolor = 'black')
ax.set_xlabel('# of classes')
ax.set_ylabel('Adjusted Mutual Information Score')
ax.set_ylim(0,1.1)
ax.set_xlim(0,10)
for i, txt in enumerate(descrip):
    ax.annotate(txt, (x[i], y[i]))
plt.savefig(resdir+'/'+data1+'/'+data1+'_AMI_number_classes_seurat_ref.pdf')
plt.clf()
