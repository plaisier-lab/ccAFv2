##########################################################
## ccAFv2:  Plotting CV classification methods          ##
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

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import scipy
from scipy import stats

##################
# Set up section #
##################

data1 = 'U5'
tags = ['SVMrej', 'RFpy', 'KNN', 'ACTINN','ccAFv2']
resdir = 'results/classification_method_comparison_redo_with_10_perc_holdout'
savedir = 'plots'
toPlot = 'f1-score'

## Load datasets and concatenate
cReport = {}
for tag1 in tags:
    if tag1 == 'ccAFv2':
        tmp1 = pd.read_csv(resdir+'/'+tag1+'_CV_classification_report_861_U5_102423_fc_25_5_layers_600_200.csv')
    else:
        tmp1 = pd.read_csv(resdir+'/'+tag1+'_CV_classification_report_861.csv')
    cReport[tag1] = tmp1

cReports = pd.concat([cReport[tags[0]], cReport[tags[1]], cReport[tags[2]], cReport[tags[3]], cReport[tags[4]]], axis = 0)
cReports.rename(columns={'Unnamed: 0':'Cell Cycle State'}, inplace=True)


########################
# Plotting f1-score   #
#######################
def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


sns.set(style="whitegrid", font_scale=3)
fig, ax = plt.subplots(figsize=(40,20))
sns.boxplot(hue="Classifier", y=toPlot, x = "Cell Cycle State", data=cReports, linewidth=1, order=['Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'])
ax.set_ylim(0,1)
ax.set(ylabel=toPlot)
plt.savefig(resdir+'/'+data1+'_'+toPlot+'_030724.pdf')
plt.clf()


# stats
stats = {}
for classifier1 in ['SVMrej', 'RFpy', 'KNN', 'ACTINN']:
    print(classifier1)
    stats[classifier1] = {}
    tmp = cReports[cReports['Classifier'] == classifier1]
    ref = cReports[cReports['Classifier'] == 'ccAFv2']
    for state1 in tmp['Cell Cycle State'].unique():
        print(state1)
        comp = tmp[tmp['Cell Cycle State'] == state1]['f1-score']
        ccafv2 = ref[ref['Cell Cycle State'] == state1]['f1-score']
        stats[classifier1][state1] = scipy.stats.ttest_ind(comp,ccafv2)[1]


df = pd.DataFrame(stats)
df.to_csv(resdir+'/ccAFv2_significantly_different_for_each_classifier_for_each_state_stats.csv')
