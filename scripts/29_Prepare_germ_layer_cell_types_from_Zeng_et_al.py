#docker run -it -v '/home/soconnor/old_home/ccNN/:/files' cplaisier/ccnn

import scanpy as sc
import pandas as pd
import numpy as np
import pickle

# read in data
bm = sc.read_h5ad('testData/GSE155121/GSE155121_human_data_raw.h5ad')
metadata = pd.read_csv('testData/GSE155121/Fig1_allweek_cluster_90_metadata.csv', index_col=0)

# subset data
adata_subset = bm[bm.obs_names.isin(metadata.index)]
adata_subset.obs = metadata
adata_subset.shape
#adata_subset.var_names_make_unique()
# (430808, 32738)

# add tsne/umap observation into obsm
adata_subset.obsm['tsne'] = np.array(adata_subset.obs[['TSNE_1', 'TSNE_2']].loc[adata_subset.obs_names])
adata_subset.obsm['umap'] = np.array(adata_subset.obs[['UMAP_1', 'UMAP_2']].loc[adata_subset.obs_names])

sc.pp.filter_cells(adata_subset, min_genes=200)
sc.pp.filter_genes(adata_subset, min_cells=3)
adata_subset.shape

# Ectoderm: neural crest (NC) cells, NSCs, intermediate progenitor cells (IPCs), and glutamatergic (Glu), GABAergic (GABA), dopaminergic (DA), and cholinergic (ACh) neurons
# Mesoderm: myoblasts, osteoblasts, primitive erythroid cells, ProMeso, lateral plate mesoderm, NMP
# Endodermal: cluster 90


for celltype1 in adata_subset.obs['Main_cell_type'].value_counts().index[1:]:
    print(celltype1)
    adata_subset[adata_subset.obs['Main_cell_type']==celltype1].write_h5ad('testData/GSE155121/new/'+celltype1+'.h5ad')
