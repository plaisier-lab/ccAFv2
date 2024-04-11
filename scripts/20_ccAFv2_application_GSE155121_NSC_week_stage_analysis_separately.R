##########################################################
## ccAFv2:  Application to GSE155121                    ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccafv2_extra

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
library(reticulate)
use_python('/usr/bin/python3')

# Install ccAFv2
devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)
mgenes = read.csv(system.file("extdata", "ccAFv2_genes.csv", package = "ccAFv2"), header = TRUE, row.names = 1)[, paste0('human_ensembl')]

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
#ccAFv2_order = c('M/Early G1', 'G2/M', 'S/G2', 'S', 'Late G1', 'G1', 'Neural G0')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")

# Set up directory structure
#setwd("files/")
resdir = 'data/GSE155121'
tag1 = 'NSC'
resdir2 = file.path(resdir, tag1)

# Load data
seurat1 <- LoadH5Seurat(file.path(resdir2,'nsc.h5seurat'))
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")
nCount = colSums(x = seurat1, slot = "counts")  # nCount_RNA
seurat1 <- AddMetaData(seurat1, nCount, col.name = "nCount_RNA")

scQC_and_norm = function(data, week_stage, v1 = 2000, v2 = 40000, h1 = 0.001, h2 = 0.1, save_dir = resdir2){
  seurat1 = data
  ws1 = week_stage
  cat('\n',week_stage,'\n')
  # Create new folders to save results
  resdir3 = file.path(save_dir, ws1)
  resdir4 = file.path(resdir3, 'analysis_output')
  resdir5 = file.path(resdir3, 'seurat_objects')
  dir.create(file.path(resdir3), showWarnings = FALSE)
  dir.create(file.path(resdir4), showWarnings = FALSE)
  dir.create(file.path(resdir5), showWarnings = FALSE)
  seurat_subset = subset(seurat1, subset = week_stage == ws1)
  pdf(file.path(resdir4, paste0(ws1, '_QC_plot_to_choose_cutoffs.pdf')))
  plot(seurat_subset@meta.data$nCount_RNA, seurat_subset@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  dev.off()
  cat('\n', dim(seurat_subset),'\n')
  # Quality control filtering
  keep.detect = which(seurat_subset@meta.data$percent.mito < h2 & seurat_subset@meta.data$percent.mito > h1 & seurat_subset@meta.data$nCount_RNA < v2 & seurat_subset@meta.data$nCount_RNA > v1)
  seurat_subset = subset(seurat_subset, cells=colnames(seurat_subset)[keep.detect])
  cat('\n', dim(seurat_subset),'\n')
  cat('Normalization \n')
  seurat_subset2 = SCTransform(seurat_subset, verbose = FALSE)
  seurat_subset2 = CellCycleScoring(object=seurat_subset2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
  phase_calls = seurat_subset2$Phase
  # Change to Ensembl ID
  seurat_subset@assays$RNA@data@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
  seurat_subset@assays$RNA@counts@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
  seurat_subset$Phase = phase_calls
  seurat_subset = SCTransform(seurat_subset, verbose = FALSE, return.only.var.genes = FALSE)
  seurat_subset = PredictCellCycle(seurat_subset, do_sctransform=FALSE)
  seurat_subset$ccAFv2 = as.character(seurat_subset$ccAFv2)
  seurat_subset$Main_cell_type = as.character(seurat_subset$Main_cell_type)
  seurat_subset$week_stage = as.character(seurat_subset$week_stage)
  saveRDS(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.rds')))
  #SaveLoom(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.loom')), verbose = FALSE, overwrite = TRUE)
  SaveH5Seurat(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.h5Seurat')), overwrite = TRUE)
  Convert(file.path(resdir5, paste0(ws1, '_normalized_ensembl.h5Seurat')), dest = 'h5ad', overwrite = TRUE)
}

scQC_and_norm(data = seurat1, week_stage = 'W3-1', h2 = 0.04, v1 = 2000)
scQC_and_norm(data = seurat1, week_stage = 'W4-1', h2 = 0.08, v1 = 4000)
scQC_and_norm(data = seurat1, week_stage = 'W4-2', h2 = 0.06, v1 = 4000, v2 = 32000)
scQC_and_norm(data = seurat1, week_stage = 'W4-3', h2 = 0.06, v1 = 4000, v2 = 28000)
scQC_and_norm(data = seurat1, week_stage = 'W5-1', h2 = 0.05, v1 = 4000, v2 = 32000)
scQC_and_norm(data = seurat1, week_stage = 'W5-2', h1 = 0.005, h2 = 0.14, v1 = 4000, v2 = 30000)
scQC_and_norm(data = seurat1, week_stage = 'W5-3', h1 = 0.005, h2 = 0.06, v1 = 4000, v2 = 35000)
scQC_and_norm(data = seurat1, week_stage = 'W6-1', h2 = 0.04, v1 = 4000, v2 = 32000)
scQC_and_norm(data = seurat1, week_stage = 'W7-1', h2 = 0.05, v1 = 4000, v2 = 30000)
scQC_and_norm(data = seurat1, week_stage = 'W8-1', h2 = 0.04, v1 = 4000, v2 = 23000)
scQC_and_norm(data = seurat1, week_stage = 'W9-1', h2 = 0.04, v1 = 4000, v2 = 23000)
scQC_and_norm(data = seurat1, week_stage = 'W9-2', h2 = 0.05, v1 = 4000, v2 = 25000)
scQC_and_norm(data = seurat1, week_stage = 'W12-1', h2 = 0.04, v1 = 4000, v2 = 20000)


datas_ccSeurat = c()
datas_ccAFv2 = c()
for (ws1 in c('W3-1', 'W4-1', 'W4-2', 'W4-3', 'W5-1', 'W5-2', 'W5-3', 'W6-1', 'W7-1', 'W8-1', 'W9-1', 'W9-2', 'W12-1')){
  resdir3 = file.path(resdir2, ws1)
  resdir4 = file.path(resdir3, 'analysis_output')
  resdir5 = file.path(resdir3, 'seurat_objects')
  tmp = readRDS(file.path(resdir5, paste0(ws1, '_normalized_ensembl.rds')))
  # Order  ccSeurat calls
  sub1 = ccSeurat_order %in% factor(tmp$Phase)
  tmp$Phase <- factor(tmp$Phase, levels = ccSeurat_order[sub1])
  df1 = data.frame(table(tmp$Phase))
  rownames(df1) = df1$Var1
  df2 = df1['Freq']/dim(tmp)[2]
  df2$Phase = rownames(df2)
  datas_ccSeurat[[ws1]] = df2
  # Order ccAFv2 calls
  sub2 = ccAFv2_order %in% factor(tmp$ccAFv2)
  tmp$ccAFv2 <- factor(tmp$ccAFv2, levels = ccAFv2_order[sub2])
  df3 = data.frame(table(tmp$ccAFv2))
  rownames(df3) = df3$Var1
  df4 = df3['Freq']/dim(tmp)[2]
  df4$ccAFv2 = rownames(df4)
  datas_ccAFv2[[ws1]] = df4
}

# ccSeurat pie plots
melted = melt(datas_ccSeurat)
melted$Phase = factor(melted$Phase, levels = ccSeurat_order[sub1])
pdf(file.path(resdir2, 'ccSeurat_pie_charts.pdf'), width = 8, height = 6)
ggplot(data=melted, aes(x=" ", y=value, group=Phase, fill=Phase)) + scale_fill_manual(values=ccSeurat_colors) +
         geom_bar(width = 1, stat = "identity") +
         coord_polar("y", start=0) +
         facet_grid(.~ L1) +theme_void() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.position='bottom') +
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
dev.off()

# ccAFv2 pie plots
melted = melt(datas_ccAFv2)
melted$ccAFv2 = factor(melted$ccAFv2, levels = ccAFv2_order[sub1])
pdf(file.path(resdir2, 'ccAFv2_pie_charts.pdf'), width = 8, height = 6)
ggplot(data=melted, aes(x=" ", y=value, group=ccAFv2, fill=ccAFv2)) + scale_fill_manual(values=ccAFv2_colors[sub1]) +
         geom_bar(width = 1, stat = "identity") +
         coord_polar("y", start=0) +
         facet_grid(.~ L1) +theme_void() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.position='bottom') +
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
dev.off()
