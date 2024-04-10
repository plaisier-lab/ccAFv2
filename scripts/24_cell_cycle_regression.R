#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccafv2_extra

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(SeuratWrappers)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(aricode)
library(tricycle)
library(peco)
library(doParallel)
library(foreach)
library(mclust)
library(cluster)
library(TSP)
library(scran)
library(yarrr)
use_python('/usr/bin/python3')

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)

install.packages('infotheo')
library('infotheo')

ccAF_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
ccAF_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")

#------------------------------------------------------
# Read in data section / set up seurat object / QC
#---------------------------------------------------

# Set working directory
setwd("files/")
mgenes = read.csv("model/final/ccAFv2_genes_full_dataset_102423_fc_25_v2_layers_600_200.csv")[,2]

# Function for variance explained
var_exp = function(data) {
    data.pca = princomp(data)
    return((data.pca$sdev^2)[1]/sum(data.pca$sdev^2))
}

# Load in marker genes per cluster
mgenes_per_cluster = read.csv('mgenes_by_cluster.csv', row.names = 'X')
cluster_mgenes = list()
num_mgenes = list()
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  mgenes_cluster = mgenes_per_cluster[mgenes_per_cluster$cluster == cluster1,]['human_ensembl']
  cluster_mgenes[[cluster1]] = mgenes_cluster$human_ensembl
  mgenes_cluster_num = length(mgenes_cluster$human_ensembl)
  # Tabulate number of marker genes per cluster
  num_mgenes[[cluster1]] = mgenes_cluster_num
}

#-------------------
# Controls
#-------------------
#-----------
# ccAFv2
#-----------
# Load data and set default assay as 'RNA'
data1 = readRDS('data/normalized/final/U5_normalized_ensembl.rds')
DefaultAssay(data1) = 'RNA'
data1 = SCTransform(data1, return.only.var.genes=FALSE)
# Apply ccAFv2
data1 = PredictCellCycle(data1, do_sctransform=FALSE)
sub1 = ccAFv2_order %in% factor(data1$ccAFv2)
data1$ccAFv2 <- factor(data1$ccAFv2, levels = ccAFv2_order[sub1])

# Add in cluster marker genes add module scores
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  cat('\n',cluster1,'\n')
  data1 = AddModuleScore(data1, features = list(cluster_mgenes[[cluster1]]), name = file.path(paste0(cluster1, '_expr')))
}

# Control
con_var_exp = list()
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  data2 = data1
  data_sub = data2[cluster_mgenes[[cluster1]],]
  data.pca = princomp(t(data_sub[['SCT']]@scale.data))
  con_var_exp[[cluster1]] = (data.pca$sdev^2)[1]/sum(data.pca$sdev^2)
}

# Regress out ccAFv2
data2 = data1
DefaultAssay(data2) = 'RNA'
# Module Scores
#data2 <- SCTransform(data2, vars.to.regress = c("Neural.G0_expr1", "G1_expr1", "S_expr1", "Late.G1_expr1", "M.Early.G1_expr1","S.G2_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)
#data2 <- SCTransform(data2, vars.to.regress = c("S_expr1", "Late.G1_expr1", "M.Early.G1_expr1","S.G2_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)
data2 <- SCTransform(data2, vars.to.regress = c("S_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)

# Control regressed
regress_var_exp = list()
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  data3 = data2
  data_sub = data3[cluster_mgenes[[cluster1]],]
  data.pca = princomp(t(data_sub[['SCT']]@scale.data))
  regress_var_exp[[cluster1]] = (data.pca$sdev^2)[1]/sum(data.pca$sdev^2)
}

compare1 = cbind(t(data.frame(con_var_exp)), t(data.frame(regress_var_exp)))
colnames(compare1) = c('con', 'ccAFv2_reg')


#-------------
# ccSeurat
#-------------
# Data preparation
data1 = readRDS('data/normalized/final/U5_normalized_ensembl.rds')
DefaultAssay(data1) = 'RNA'
data1 = SCTransform(data1, return.only.var.genes=FALSE)
# Read in ccSeurat info
s.score = read.csv('U5_ccSeurat_S_Score.csv', row.names = 'X')
g2m.score = read.csv('U5_ccSeurat_G2M_Score.csv', row.names = 'X')
ccSeurat_call = read.csv('U5_ccSeurat_Phase.csv', row.names = 'X')
# Put into object metadata
data1$S.Score = s.score$x
data1$G2M.Score = g2m.score$x
data1$Phase = ccSeurat_call$x
# Organize calls
data1$Phase <- factor(data1$Phase, levels = ccSeurat_order)

# Regress out ccSeurat
DefaultAssay(data1) = 'RNA'
data2 <- SCTransform(data1, vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes=FALSE)
seurat_regress_var_exp = list()
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  data3 = data2
  data_sub = data3[cluster_mgenes[[cluster1]],]
  data.pca = princomp(t(data_sub[['SCT']]@scale.data))
  seurat_regress_var_exp[[cluster1]] = (data.pca$sdev^2)[1]/sum(data.pca$sdev^2)
}

compare2 = cbind(compare1, t(data.frame(seurat_regress_var_exp)))
colnames(compare2) = c('con', 'ccAFv2_reg', 'ccSeurat_reg')
rownames(compare2) = c('G1', 'S', 'Neural G0', 'Late G1', 'M/Early G1', 'S/G2', 'G2/M')
write.csv(data.frame(compare2), 'ccAFv2_marker_genes_variance_explained_ccAFv2_and_ccseurat_addmodulescore.csv')


#---------------
# Permutations
#---------------

#-----------
# ccAFv2
#-----------

 Load data and set default assay as 'RNA'
data1 = readRDS('data/normalized/final/U5_normalized_ensembl.rds')
DefaultAssay(data1) = 'RNA'
data1 = SCTransform(data1, return.only.var.genes=FALSE)
# Apply ccAFv2
data1 = PredictCellCycle(data1, do_sctransform=FALSE)
sub1 = ccAFv2_order %in% factor(data1$ccAFv2)
data1$ccAFv2 <- factor(data1$ccAFv2, levels = ccAFv2_order[sub1])

# Subset genes to union of top 3000 variable genes and 861 ccAFv2 genes
genes_sub = Reduce(union, list(mgenes, data1@assays$SCT@var.features))
data1 = data1[genes_sub,]

for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  cat('\n',cluster1,'\n')
  data1 = AddModuleScore(data1, features = list(cluster_mgenes[[cluster1]]), name = file.path(paste0(cluster1, '_expr')))
}

# Regress out ccAFv2
data2 = data1
DefaultAssay(data2) = 'RNA'
#data2 <- SCTransform(data2, vars.to.regress = c("S_expr1", "Late.G1_expr1", "M.Early.G1_expr1","S.G2_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)
#data2 <- SCTransform(data2, vars.to.regress = c("S_expr1", "Late.G1_expr1", "M.Early.G1_expr1","S.G2_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)
data2 <- SCTransform(data2, vars.to.regress = c("S_expr1", "G2.M_expr1"), return.only.var.genes=FALSE)

# Run for regressed ccAFv2:
mean_ccAFv2_regressed_permuted_pvals = list()
std_ccAFv2_regressed_permuted_pvals = list()
ccAFv2_regressed_permuted_pvalue = list()
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  cat('\n',cluster1,'\n')
  data3 = data2
  numSamples = num_mgenes[[cluster1]]
  allInds = 1:length(rownames(data3))
  sampWOreplace = replicate(nfolds,sample(allInds,numSamples,replace = FALSE))
  ccAFv2_regressed_permuted_pvalue[[cluster1]] = 0
  permuted_var_exp = sapply(1:nfolds, function(k) { var_exp(t((data3[sampWOreplace[,k],][['SCT']]@scale.data))) } )
  mean_ccAFv2_regressed_permuted_pvals[[cluster1]] = mean(permuted_var_exp)
  std_ccAFv2_regressed_permuted_pvals[[cluster1]] = sd(permuted_var_exp)
  ccAFv2_regressed_permuted_pvalue[[cluster1]] = sum(permuted_var_exp >= data.frame(compare2)[cluster1,]$ccAFv2_reg)/nfolds
}

cf = t(data.frame(ccAFv2_regressed_permuted_pvalue))


#-----------
# ccSeurat
#-----------
# ccSeurat
# Load in data
data1 = readRDS('data/normalized/final/U5_normalized_ensembl.rds')
DefaultAssay(data1) = 'RNA'
data1 = SCTransform(data1, return.only.var.genes=FALSE)
# Read in ccSeurat info
s.score = read.csv('U5_ccSeurat_S_Score.csv', row.names = 'X')
g2m.score = read.csv('U5_ccSeurat_G2M_Score.csv', row.names = 'X')
ccSeurat_call = read.csv('U5_ccSeurat_Phase.csv', row.names = 'X')
# Put into object metadata
data1$S.Score = s.score$x
data1$G2M.Score = g2m.score$x
data1$Phase = ccSeurat_call$x
# Organize calls
data1$Phase <- factor(data1$Phase, levels = ccSeurat_order)

# Regress out ccSeurat
DefaultAssay(data1) = 'RNA'
data2 <- SCTransform(data1, vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes=FALSE)

# Subset genes to union of top 3000 variable genes and 861 ccAFv2 genes
genes_sub = Reduce(union, list(mgenes, data2@assays$SCT@var.features))
data2 = data2[genes_sub,]

# Run for regressed ccSeurat:
mean_seurat_regressed_permuted_pvals = list()
std_seurat_regressed_permuted_pvals = list()
seurat_regressed_permuted_pvalue = list()
nfolds = 1000
for (cluster1 in unique(mgenes_per_cluster['cluster'])$cluster){
  cat('\n',cluster1,'\n')
  data3 = data2
  numSamples = num_mgenes[[cluster1]]
  allInds = 1:length(rownames(data3))
  sampWOreplace = replicate(nfolds,sample(allInds,numSamples,replace = FALSE))
  seurat_regressed_permuted_pvalue[[cluster1]] = 0
  permuted_var_exp = sapply(1:nfolds, function(k) { var_exp(t((data3[sampWOreplace[,k],][['SCT']]@scale.data))) } )
  mean_seurat_regressed_permuted_pvals[[cluster1]] = mean(permuted_var_exp)
  std_seurat_regressed_permuted_pvals[[cluster1]] = sd(permuted_var_exp)
  seurat_regressed_permuted_pvalue[[cluster1]] = sum(permuted_var_exp >= data.frame(compare2)[cluster1,]$ccSeurat_reg)/nfolds
}

df = cbind(cf, seurat_regressed_permuted_pvalue)
colnames(df) = c('ccAFv2_reg', 'seurat_reg')

# Save out empirical pvalues
#write.csv(df, 'permutated_and_regressed_ccAFv2_pvalues_3365_genes_1000_iters_addmodulescore_all_cycling_regressed.csv')
write.csv(df, 'permutated_and_regressed_ccAFv2_pvalues_3365_genes_1000_iters_addmodulescore_s_and_g2m_regressed.csv')

# Save out mean and standard dev
cf2 = cbind(t(data.frame(mean_ccAFv2_regressed_permuted_pvals)), t(data.frame(std_ccAFv2_regressed_permuted_pvals)))
colnames(cf2) = c('permutated_ccAFv2_reg_mean', 'permutated_ccAFv2_reg_stdev')
#write.csv(cf2, 'permutated_and_regressed_ccAFv2_pvalues_mean_and_stdev_3365_genes_1000_iters_addmodulescore_all_cycling.csv')
write.csv(cf2, 'permutated_and_regressed_ccAFv2_pvalues_mean_and_stdev_3365_genes_1000_iters_addmodulescore_only_s_g2m.csv')
