##########################################################
## ccAFv2:  Plot Zeng et al. cell types                 ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/testData/GSE155121/new:/files' cplaisier/ccnn

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import functools
import os
import itertools

#tags = ['Notochord', 'Endoderm', 'IPC', 'Endothelium', 'Blood', 'NMP', 'Erythroid', 'Immune Cell', 'Cholin N', 'Epithelium', 'LPM', 'Neural Crest', 'Osteoblast', 'Chondroblast', 'Glu N', 'Myoblast', 'DA N', 'GABA N', 'MesoPro']
tags = ['Endoderm', 'IPC', 'Blood', 'NMP', 'Erythroid', 'Cholin N', 'LPM', 'Neural Crest', 'Osteoblast', 'Chondroblast', 'Glu N', 'Myoblast', 'DA N', 'GABA N', 'MesoPro']
week_stages = ['W3-1', 'W4-1', 'W4-2', 'W4-3', 'W5-1', 'W5-2', 'W5-3', 'W6-1', 'W7-1', 'W8-1', 'W9-1', 'W9-2', 'W12-1']
ccAFv2_colors = ["#d9a428", "#f37f73", "#1fb1a9", "#8571b2", "#db7092", "#3db270" , "#6d90ca",  "#D3D3D3"]

# Read in data
allData = pd.DataFrame()
lay = 'NA'
for tag in tags:
    if tag in ['Endoderm']:
        lay = 'Endoderm'
    if tag in ['Myoblast', 'Osteoblast', 'Chondroblast', 'Erythroid', 'LPM', 'NMP', 'Blood', 'MesoPro']:
        lay = 'Mesoderm'
    if tag in ['DA N', 'GABA N', 'Neural Crest', 'Cholin N', 'IPC', 'Glu N']:
        lay = 'Ectoderm'
    for ws in week_stages:
        if os.path.isfile(tag+'/'+ws+'/analysis_output/'+tag+'_ccAFv2_calls.csv'):
            data1 = pd.read_csv(tag+'/'+ws+'/analysis_output/'+tag+'_ccAFv2_calls.csv', index_col = 0)
            data1['Main_cell_type'] = tag
            data1['Week_stage'] = ws
            data1['Dermal_layer'] = lay
            allData = allData.append(data1)


#allData.to_csv('ccAFv2_calls_for_Zeng.csv')
allData = pd.read_csv('ccAFv2_calls_for_Zeng.csv', index_col = 0)

# Plot
sns.histplot(
    data=allData,
    x="Dermal_layer", hue="x",
    multiple="fill", stat="proportion",
    discrete=True, shrink=.8,
    hue_order = ['Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown'],
    palette = {'Neural G0': '#d9a428', 'G1': '#f37f73', 'Late G1':'#1fb1a9', 'S': '#8571b2', 'S/G2':'#db7092', 'G2/M': '#3db270', 'M/Early G1': '#6d90ca', 'Unknown': '#D3D3D3'}
)
plt.xticks(rotation=90)
plt.savefig('Zeng_ccAFv2_dermal.pdf')


# Find common Neural G0 endodermal marker genes
endo_mgenes = {}
sep_tags = ['Endoderm']
for tag1 in sep_tags:
    endo_mgenes[tag1] = []
    for ws1 in week_stages:
        if os.path.isfile(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv'):
            data1 = pd.read_csv(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv', index_col = 0)
            data1 = data1[data1['cluster'] == 'Neural G0']
            data2 = list(data1[data1['p_val_adj']<=0.05].index)
            endo_mgenes[tag1].append(data2)

endoderm_mgenes = set(list(itertools.chain(*endo_mgenes['Endoderm'])))
pd.DataFrame(list(endoderm_mgenes)).to_csv('endoderm_mgenes.csv')
df = pd.DataFrame(list(endoderm_mgenes))
endoderm_common_mgenes = df

# Find common Neural G0 ectodermal marker genes
ecto_mgenes = {}
sep_tags = ['DA N', 'GABA N', 'Neural Crest', 'Cholin N', 'IPC', 'Glu N']
for tag1 in sep_tags:
    ecto_mgenes[tag1] = []
    for ws1 in week_stages:
        if os.path.isfile(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv'):
            data1 = pd.read_csv(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv', index_col = 0)
            data1 = data1[data1['cluster'] == 'Neural G0']
            data2 = list(data1[data1['p_val_adj']<=0.05].index)
            ecto_mgenes[tag1].append(data2)

da_n_mgenes = set(list(itertools.chain(*ecto_mgenes['DA N'])))
pd.DataFrame(list(da_n_mgenes)).to_csv('da_n_mgenes.csv')
glu_n_mgenes = set(list(itertools.chain(*ecto_mgenes['Glu N'])))
pd.DataFrame(list(glu_n_mgenes)).to_csv('glu_n_mgenes.csv')
neural_crest_mgenes = set(list(itertools.chain(*ecto_mgenes['Neural Crest'])))
pd.DataFrame(list(neural_crest_mgenes)).to_csv('neural_crest_mgenes.csv')
cholin_n_mgenes = set(list(itertools.chain(*ecto_mgenes['Cholin N'])))
pd.DataFrame(list(cholin_n_mgenes)).to_csv('cholin_n_mgenes.csv')
ipc_mgenes = set(list(itertools.chain(*ecto_mgenes['IPC'])))
pd.DataFrame(list(ipc_mgenes)).to_csv('ipc_mgenes.csv')
gaba_n_mgenes = set(list(itertools.chain(*ecto_mgenes['GABA N'])))
pd.DataFrame(list(gaba_n_mgenes)).to_csv('gaba_n_mgenes.csv')

tmp = da_n_mgenes & glu_n_mgenes & neural_crest_mgenes & cholin_n_mgenes & ipc_mgenes & gaba_n_mgenes
#ENSG00000046653 in all

tmp2 = list(da_n_mgenes) + list(glu_n_mgenes) + list(neural_crest_mgenes) + list(cholin_n_mgenes) + list(ipc_mgenes) + list(gaba_n_mgenes)
df = pd.DataFrame(tmp2)
df2 = pd.DataFrame(df.value_counts())
df2.to_csv('ectodermal_common_mgenes.csv')
ectodermal_common_mgenes = pd.read_csv('ectodermal_common_mgenes.csv', index_col =0)
df3 = df2[df2[0]>=5]

# Find common Neural G0 mesodermal marker genes
meso_mgenes = {}
sep_tags = ['Myoblast', 'Osteoblast', 'Chondroblast', 'Erythroid', 'LPM', 'NMP', 'Blood', 'MesoPro']
for tag1 in sep_tags:
    meso_mgenes[tag1] = []
    for ws1 in week_stages:
        if os.path.isfile(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv'):
            data1 = pd.read_csv(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv', index_col = 0)
            data1 = data1[data1['cluster'] == 'Neural G0']
            data2 = list(data1[data1['p_val_adj']<=0.05].index)
            meso_mgenes[tag1].append(data2)

myoblast_mgenes = set(list(itertools.chain(*meso_mgenes['Myoblast'])))
pd.DataFrame(list(myoblast_mgenes)).to_csv('myoblast_mgenes.csv')
chondroblast_mgenes = set(list(itertools.chain(*meso_mgenes['Chondroblast'])))
pd.DataFrame(list(chondroblast_mgenes)).to_csv('chondroblast_mgenes.csv')
osteoblast_mgenes = set(list(itertools.chain(*meso_mgenes['Osteoblast'])))
pd.DataFrame(list(osteoblast_mgenes)).to_csv('osteoblast_mgenes.csv')
lpm_mgenes = set(list(itertools.chain(*meso_mgenes['LPM'])))
pd.DataFrame(list(lpm_mgenes)).to_csv('lpm_mgenes.csv')
nmp_mgenes = set(list(itertools.chain(*meso_mgenes['NMP'])))
pd.DataFrame(list(nmp_mgenes)).to_csv('nmp_mgenes.csv')
erythroid_mgenes = set(list(itertools.chain(*meso_mgenes['Erythroid'])))
pd.DataFrame(list(erythroid_mgenes)).to_csv('erythroid_mgenes.csv')
blood_mgenes = set(list(itertools.chain(*meso_mgenes['Blood'])))
pd.DataFrame(list(blood_mgenes)).to_csv('blood_mgenes.csv')
mesopro_mgenes = set(list(itertools.chain(*meso_mgenes['MesoPro'])))
pd.DataFrame(list(mesopro_mgenes)).to_csv('mesopro_mgenes.csv')


myoblast_mgenes = set(list(pd.read_csv('myoblast_mgenes.csv', index_col = 0)['0']))
chondroblast_mgenes = set(list(pd.read_csv('chondroblast_mgenes.csv', index_col = 0)['0']))
osteoblast_mgenes  = set(list(pd.read_csv('osteoblast_mgenes.csv', index_col = 0)['0']))
lpm_mgenes  = set(list(pd.read_csv('lpm_mgenes.csv', index_col = 0)['0']))
nmp_mgenes  = set(list(pd.read_csv('nmp_mgenes.csv', index_col = 0)['0']))
erythroid_mgenes  = set(list(pd.read_csv('erythroid_mgenes.csv', index_col = 0)['0']))
blood_mgenes  = set(list(pd.read_csv('blood_mgenes.csv', index_col = 0)['0']))
mesopro_mgenes  = set(list(pd.read_csv('mesopro_mgenes.csv', index_col = 0)['0']))

tmp = myoblast_mgenes & chondroblast_mgenes & osteoblast_mgenes & lpm_mgenes & blood_mgenes & mesopro_mgenes & erythroid_mgenes & nmp_mgenes

tmp2 = list(myoblast_mgenes) + list(chondroblast_mgenes) + list(osteoblast_mgenes) + list(lpm_mgenes) + list(blood_mgenes) + list(mesopro_mgenes) + list(erythroid_mgenes) + list(nmp_mgenes)
df = pd.DataFrame(tmp2)
df2 = pd.DataFrame(df.value_counts())
df2.to_csv('mesodermal_common_mgenes.csv')

mesodermal_common_mgenes = pd.read_csv('mesodermal_common_mgenes.csv', index_col = 0)
endoderm_common_mgenes = pd.read_csv('endoderm_mgenes.csv', index_col = 0)
ectodermal_common_mgenes = pd.read_csv('ectodermal_common_mgenes.csv', index_col = 0)
m1 = list(mesodermal_common_mgenes[mesodermal_common_mgenes['0']>=4].index)
en = list(endoderm_common_mgenes['0'])
ec = list(ectodermal_common_mgenes[ectodermal_common_mgenes['0']>=3].index)


all1 = m1 + en + ec
all_df = pd.DataFrame(all1)
all_df.value_counts()

tmp = pd.DataFrame(all_df.value_counts())
tmp2 = tmp.rename(columns = {0:'count'})
tmp2.reset_index(inplace=True)
tmp3 = tmp2[tmp2['count']>=2]
common_genes = list(tmp3[0])




# common Neural G0 marker genes
endo_mgenes = {}
meso_mgenes = {}
ecto_mgenes = {}

tags = ['GABA N']
lays = ['Ectoderm']
ws = ['W4-2', 'W6-1', 'W7-1', 'W8-1', 'W9-1', 'W9-2', 'W12-1']

for tag1 in tags:
    ecto_mgenes[tag1] = pd.DataFrame()
    for ws1 in ws:
        data1 = pd.read_csv(tag1+'/'+ws1+'/analysis_output/'+ws1+'_scTransform_Markers_together.csv', index_col = 0)
        data1 = data1[data1['cluster'] == 'Neural G0']
        data2 = data1[data1['p_val_adj']<=0.05]
        ecto_mgenes[tag1] = pd.concat([ecto_mgenes[tag1], data2])

ecto_mgenes[tag1].reset_index(inplace = True)
data2 = pd.DataFrame(ecto_mgenes[tag1].groupby('index')['avg_log2FC'].max())

#common_genes_fc = pd.DataFrame(index = common_genes, columns = ['Endoderm', 'Myoblast', 'Chondroblast', 'Osteoblast', 'LPM', 'NMP', 'Erythroid', 'Blood', 'MesoPro', 'DA N', 'Glu N', 'Neural Crest', 'Cholin N', 'IPC', 'GABA N'])
for gene1 in common_genes:
    if gene1 in data2.index:
        fc1 = round(data2.loc[gene1]['avg_log2FC'],3)
        common_genes_fc.loc[gene1][tag1] = fc1
    else:
        common_genes_fc.loc[gene1][tag1] = 0

common_genes_fc.to_csv('Zeng_common_G0_genes_across_main_cell_types.csv')
common_genes_fc = pd.read_csv('Zeng_common_G0_genes_across_main_cell_types.csv')
common_genes_fc = common_genes_fc.rename(columns={'Unnamed: 0':'Marker gene'})
common_genes_fc.set_index('Marker gene', inplace=True)

plt.rcParams['figure.figsize']=(10,10)
sns.heatmap(common_genes_fc.T, linewidth=.5, cmap="Reds")
plt.savefig('Zeng_common_g0_mgenes_v3.pdf')
plt.close()


#df = pd.DataFrame(index = range(1,22), columns = {'Endo', 'Meso', 'Ecto'})
df.loc[21]['Endo'] =0
df.loc[21]['Meso'] =7
df.loc[21]['Ecto'] =5

df.reset_index(inplace=True)
df.rename(columns = {'index':'gene'}, inplace=True)

plt.bar(df['gene'], df['Ecto'], color='b')
plt.bar(df['gene'], df['Meso'], bottom = df['Ecto'], color='r')
plt.bar(df['gene'], df['Endo'], bottom = df['Ecto'], color='g')
plt.savefig('Zeng_mgenes_dermal_stacked_v2.pdf')

df[['Ecto', 'Meso', 'Endo']].plot(kind='bar', stacked=True, color=['skyblue', 'red', 'green'])
plt.savefig('Zeng_NG0_mgenes_dermal_stacked_bar.pdf')
