##########################################################
## ccAFv2: Testing ccAFv2 600_200 model                 ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccafv2_extra

#------------------------
# Set up / imports
#-----------------------

library(dplyr)
library(readr)
library(data.table)
library(grid)
library(gridExtra)
library(Seurat)
library(SeuratDisk)
library(keras)
library(aricode)
use_python('/usr/bin/python3')

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)

# Set working directory
setwd("files/")

# Set output directory
save_dir = 'testData/all_test_data_objects'

# Cell cycle state order and colors for plotting
ccSeurat_order = c('G1', 'S', 'G2M')
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")
ccAFv2_order = c('G1/other','Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("G1/other"= "#9aca3c","Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

#-----------------------------------------------
# Load data and calculate AMI/Cells Predicted
#----------------------------------------------

# Initialize AMI and Cells Predicted
AMI1 = list()
CellsPredict = list()

# Data to use for ccAFv2 model testing
datas = c('U5', 'LGG275_GF', 'BT322', 'BT324', 'BT326', 'BT333', 'BT363', 'BT368', 'BT363_tumor', 'BT368_tumor', 'GSE155121')
for(datas1 in datas){
  # Load in data
  cat('\n', datas1,'\n')
  if(datas1 == 'U5'){
    resdir = file.path('data/normalized/final')
    resdir2 = file.path('testData/U5')
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccseurat_calls_101023.csv')), row.names = 'X')
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if (datas1 == 'LGG275_GF'){
    resdir = file.path(paste0('testData/SCDataHugnotIGF/', strsplit(datas1, split = "_")[[1]][1], '/', datas1, '/seurat_objects'))
    resdir2 = file.path(paste0('testData/SCDataHugnotIGF/', strsplit(datas1, split = "_")[[1]][1], '/', datas1, '/analysis_output'))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 'X')
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if(datas1 %in% c('BT324', 'BT326', 'BT363', 'BT368')){
    resdir = file.path(paste0('testData/GSC_bam/',datas1, '/seurat_objects'))
    resdir2 = file.path(paste0('testData/GSC_bam/',datas1))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_phase.csv')), row.names = 1)
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$seurat2.Phase, col.name = 'Phase')
  }else if(datas1 %in% c('BT322','BT333')){
    resdir = file.path(paste0('testData/GSC_bam/',datas1, '/seurat_objects'))
    resdir2 = file.path(paste0('testData/GSC_bam/',datas1))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_phase.csv')), row.names = 1)
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if(datas1 %in% c('BT363_tumor', 'BT368_tumor')){
    resdir = file.path(paste0('testData/GBM_tumors/',datas1, '/seurat_objects'))
    resdir2 = file.path(paste0('testData/GBM_tumors/',datas1, '/analysis_output'))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 1)
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if(datas1 == 'W8-1'){
    resdir = file.path(paste0('testData/GSE155121/NSC/FINAL/',datas1))
    resdir2 = resdir
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_processed_010524.rds')))
  }
  # Classify with ccAFv2
  seurat2 = PredictCellCycle(seurat2)
  # Save out ccAFv2 calls
  write.csv(seurat2$ccAFv2, file.path(resdir2, paste0(datas1, '_ccAFv2_calls_R_test.csv')))
  # Save out data as RDS object
  saveRDS(seurat2, file.path(save_dir, paste0(datas1, '_normalized_with_ccAFv2_R_calls.rds')))
  # Calculate AMI
  subset1 = seurat2[,seurat2$ccAFv2 != 'Unknown']
  predlab = subset1$ccAFv2
  predlab2 = list(data.frame(predlab)$predlab)[[1]]
  truelab = subset1$Phase
  truelab2 = list(data.frame(truelab)$truelab)[[1]]
  AMI1[[datas1]] <- round(AMI(truelab2, predlab2), 2)
  CellsPredict[[datas1]] = round((dim(subset1)[2]/dim(seurat2)[2])*100, 2)
  # Prepare for plotting / order and colors
  sub1 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub1])
  sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
  # Plot
  p1 = DimPlot(seurat2, reduction = "umap", label=F, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub1])
  c2 = DimPlot(seurat2, reduction = "umap", label=F, pt.size = 1, group.by = 'ccAFv2', cols = ccAFv2_colors[sub2])
  # Save as pdf
  pdf(file.path(resdir2, paste0(datas1, '_umap_ccAFv2_phase.pdf')), width = 8, height = 10)
  lst = list(p1, c2)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA),c(2,NA)), top = "")
  dev.off()
}

# Save out AMI and Cells Predicted values for all test datasets
write.csv(data.frame(AMI1), file.path('testData/all_test_data_AMI_ccseurat_as_ref.csv'))
write.csv(data.frame(CellsPredict), file.path('testData/all_test_data_cells_predict.csv'))
