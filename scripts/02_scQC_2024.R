##########################################################
## ccAFv2: Quality control & data preparation           ##
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

# Run docker
#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccafv2_extra
#docker run -it -v '/media:/files' cplaisier/ccafv2_seurat5
#docker run -it -v '/home/soconnor/old_home:/files' cplaisier/ccafv2_seurat4

#--------------------------------
# Set up section / load packages
#--------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library("org.Hs.eg.db")
library(aricode)

# Set working directory
#setwd("files/")
resdir = 'data'

# Mitochondrial genes as ensembl IDs
mito_genes = read_csv(file.path(resdir, 'mito_genes.csv'), show_col_types = FALSE) %>% pull(mito)

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#------------------------------------------------------
# QC Function
#---------------------------------------------------

scQC = function(res_dir, tag, mito_genes, v1 = 2000, v2 = 30000, h1 = 0.001, h2 = 0.1, data_dir = 'outs/filtered_feature_bc_matrix', save_dir = 'analysis_output', obj_dir = 'seurat_objects', mt = 'MT-', symbol = F, norm_regress = F){
  cat('\n',tag,'\n')
  resdir1 = file.path(res_dir, tag)
  # Create folders
  dir.create(file.path(resdir1, save_dir), showWarnings = FALSE)
  resdir2 = file.path(resdir1, save_dir)
  dir.create(file.path(resdir1, obj_dir), showWarnings = FALSE)
  resdir3 = file.path(resdir1, obj_dir)
  #---------------------
  # Load in data
  #---------------------
  gene_column = 1
  gene_id = 'ensembl'
  if(symbol){
    gene_column = 2
    gene_id = 'gene_symbols'
  }
  data = Read10X(file.path(resdir1, data_dir), gene.column=gene_column) # column 1 is ensembl (in 10X mtx file)
  cat('Raw data: ', dim(data)[2], 'cells', dim(data)[1], 'genes \n')
  # Substitute underscores if necessary
  rownames(data) = gsub("_", "-", rownames(data))
  # Create seurat object
  seurat1 = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
  cat('Basic filter', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
  if(symbol){
    mito_genes = grep(mt, rownames(seurat1))
  }
  seurat1[['percent.mito']] = PercentageFeatureSet(seurat1, features = mito_genes)/100
  #---------------------
  # Quality control
  #---------------------
  cat('\n Quality control \n')
  # Quality control plots for choosing cutoffs
  pdf(file.path(resdir2, paste0(tag, '_QC_plot_to_choose_cutoffs.pdf')))
  plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  # Quality control filtering
  keep.detect = which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 & seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)
  seurat1 = subset(seurat1, cells=colnames(seurat1)[keep.detect])
  cat('Filtered to: ', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
  saveRDS(seurat1, file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id), '.rds')))
  # Save as new object so can go back to previous non-normalized / scaled seurat object if need too
  seurat2 = seurat1
  #------------------------------------------------------
  # Normalization with sctransform
  #---------------------------------------------------
  cat('\n Normalization \n')
  if(norm_regress){
    seurat2 = SCTransform(seurat2, verbose = FALSE, vars.to.regress = c('nCount_RNA'))
  } else {
  seurat2 = SCTransform(seurat2, verbose = FALSE)
  }
  cat('Normalized genes: ', dim(seurat2@assays$SCT@data)[1], 'features,', length(seurat2@assays$SCT@var.features), 'highly variable genes \n')
  # Classify with ccSeurat and save out as csv
  if(symbol){
    seurat2 <- CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
    write.csv(seurat2$Phase, file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')))
  }
  cat('\n Saving normalized RDS object \n')
  saveRDS(seurat2, file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id), '.rds')))
  return(seurat2)
}


#--------------------------------------
# Quality control & data preparation
#--------------------------------------

qc_ensembl = list()
qc_ensembl [['BT322']] = scQC(res_dir = 'data/GSC', tag = 'BT322', mito_genes = mito_genes, v1 = 4000, v2 = 62000, h1 = 0.009, h2 = 0.1, symbol = F, norm_regress = T)
qc_ensembl [['BT324']] = scQC(res_dir = 'data/GSC', tag = 'BT324', mito_genes = mito_genes, v1 = 5000, v2 = 40000, h1 = 0.009, h2 = 0.06, symbol = F, norm_regress = T)
qc_ensembl[['LGG275_GF']] = scQC(res_dir = 'data/LGG/LGG275', tag = 'LGG275_GF', mito_genes = mito_genes, v1 = 5000, v2 = 76000, h1 = 0.001, h2 = 0.15)
qc_ensembl[['LGG275_noGF']] = scQC(res_dir = 'data/LGG/LGG275', tag = 'LGG275_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 35000, h1 = 0.001, h2 = 0.15)
qc_ensembl[['BT237_GF']] = scQC(res_dir = 'data/LGG/BT237', tag = 'BT237_GF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.01, h2 = 0.15)
qc_ensembl[['BT237_noGF']] = scQC(res_dir = 'data/LGG/BT237', tag = 'BT237_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.17)


qc_symbol = list()
qc_symbol[['BT322']] = scQC(res_dir = 'data/GSC', tag = 'BT322', mito_genes = mito_genes, v1 = 4000, v2 = 62000, h1 = 0.009, h2 = 0.1, symbol = T, norm_regress = T)
qc_symbol[['BT324']] = scQC(res_dir = 'data/GSC', tag = 'BT324', mito_genes = mito_genes, v1 = 5000, v2 = 40000, h1 = 0.009, h2 = 0.06, symbol = T, norm_regress = T)
qc_symbol[['LGG275_GF']] = scQC(res_dir = 'data/LGG/LGG275', tag = 'LGG275_GF', mito_genes = mito_genes, v1 = 5000, v2 = 76000, h1 = 0.001, h2 = 0.15, symbol = T)
qc_symbol[['LGG275_noGF']] = scQC(res_dir = 'data/LGG/LGG275', tag = 'LGG275_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 35000, h1 = 0.001, h2 = 0.15, symbol = T)
qc_symbol[['BT237_GF']] = scQC(res_dir = 'data/LGG/BT237', tag = 'BT237_GF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.01, h2 = 0.15, symbol = T)
qc_symbol[['BT237_noGF']] = scQC(res_dir = 'data/LGG/BT237', tag = 'BT237_noGF', mito_genes = mito_genes, v1 = 6000, v2 = 90000, h1 = 0.03, h2 = 0.17, symbol = T)
