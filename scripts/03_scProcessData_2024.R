##########################################################
## ccAFv2: scRNA-seq downstream analysis                ##
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
#library(reticulate)
use_python('/usr/bin/python3')

# Set working directory
#setwd("files/")

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)
mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0('human_ensembl')]

# Some features to investigate
features1 = c('S100B', 'SOX2', 'SOX4', 'MKI67', 'APOE', 'VIM', 'CLU', 'FABP7','OLIG1','OLIG2', 'DLL3', 'HES6')
# convert to ensembl IDs
ensembl_features1 = mapIds(org.Hs.eg.db, keys = features1, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features1 = na.omit(data.frame(ensembl_features1))
ensembl_features1_plot = ensembl_features1$ensembl_features1

features2 <- c("MBP", "PLP1", "ETNPPL", "CD14","CX3CR1","PTPRC", "RBFOX3")
# convert to ensembl IDs
ensembl_features2 = mapIds(org.Hs.eg.db, keys = features2, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features2 = na.omit(data.frame(ensembl_features2))
ensembl_features2_plot = ensembl_features2$ensembl_features2

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

# Process data individually
scProcessData = function(res_dir, tag, cutoff = 0.5, assay = 'SCT', do_sctransform = TRUE, species= 'human', resolution = 0.8, save_dir = 'analysis_output', obj_dir = 'seurat_objects', clusters_to_remove = F, symbol = F, norm_regress = F){
  cat('\n',tag,'\n')
  # Set up folders
  resdir1 = file.path(res_dir, tag)
  resdir2 = file.path(resdir1, save_dir)
  resdir3 = file.path(resdir1, obj_dir)
  #---------------------------
  # Load filtered / normalized data
  #---------------------------
  gene_id = 'ensembl'
  #seurat2 = readRDS(file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id),'.rds')))
  seurat2 = readRDS(file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id),'.rds')))
  cat('\n', dim(seurat2), '\n')
  # Load ccSeurat calls
  ccseurat_calls = read.csv(file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')), row.names = 'X')
  seurat2 <- AddMetaData(seurat2, ccseurat_calls, col.name = 'Phase')
  # Order ccSeurat calls
  sub1 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub1])
  #---------------------------
  # Classify with ccAFv2
  #---------------------------
  seurat2 = PredictCellCycle(seurat2, cutoff = cutoff, assay = assay, do_sctransform = do_sctransform, species = species, gene_id = gene_id)
  # Order ccAFv2 calls
  sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
  tmp = data.frame(table(seurat2$ccAFv2))
  rownames(tmp) = tmp$Var1
  write.csv(seurat2$ccAFv2, file.path(resdir2, paste0(tag, '_ccAFv2_calls.csv')))
  write.csv(data.frame((tmp['Freq']/dim(seurat2)[2])*100), file.path(resdir2, paste0(tag, '_ccAFv2_call_frequency.csv')))
  #---------------------------
  # Normalize
  #---------------------------
  #seurat2 = SCTransform(seurat2, verbose = FALSE) # started with normalized object; do not need to redo
  # Add marker gene counts for each cell
  seurat_subset = seurat2[mgenes]
  mgene_counts = GetAssayData(seurat_subset, slot = 'counts')
  non_zero_mgenes = colSums(mgene_counts > 0)
  cat('Cells that have non-zero ccAFv2 genes: ', length(non_zero_mgenes), '\n')
  # Add as meta data column
  seurat2 = AddMetaData(seurat2, non_zero_mgenes, col.name = 'ccAFv2_mgene_counts')
  seurat2 = RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 = RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  # Check expression of non-neuronal markers
  pdf(file.path(resdir2, paste0(tag, '_check_for_non_neuronal_features.pdf')))
  print(DimPlot(seurat2, reduction = 'umap'))
  print(FeaturePlot(seurat2, reduction = 'umap', features = ensembl_features2_plot))
  dev.off()
  if(is.list(clusters_to_remove)){
    # Subset data to remove non-tumor cells
    remove1 = unlist(clusters_to_remove)
    seurat1 = seurat2
    DefaultAssay(seurat1) = 'RNA'
    Idents(object = seurat1) = "seurat_clusters"
    subset1 = subset(seurat1, idents = remove1, invert = TRUE)
    subset2 = subset1
    if(norm_regress){
      subset2 = SCTransform(subset2, verbose = FALSE, vars.to.regress = c('nCount_RNA'))
    } else {
      subset2 = SCTransform(subset2, verbose = FALSE)
    }
    subset2 = PredictCellCycle(subset2, cutoff = cutoff, assay = assay, do_sctransform = do_sctransform, species = species, gene_id = gene_id)
    subset2 = RunPCA(subset2, dims=1:30)
    subset2 = FindNeighbors(subset2, dims=1:30)
    subset2 = FindClusters(subset2, resolution = resolution)
    subset2 = RunUMAP(subset2, dims=1:30)
    seurat2 = subset2
  }
  Idents(seurat2) = seurat2$ccAFv2
  # Find cluster marker genes for each ccAFv2 class
  cluster_markers = FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE)
  cluster_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
  cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
  cluster_markers$gene = cluster_markers_genes
  write.csv(cluster_markers, file.path(resdir2, paste0(tag,'_scTransform_Markers_together.csv')))
  top10_symbol = mapIds(org.Hs.eg.db, keys = top10$gene, keytype = 'ENSEMBL', column='SYMBOL', multiVals='first')
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('Plotting UMAPs and ccAFv2 barplot \n')
  d1 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'seurat_clusters') + ggtitle('seurat_clusters')
  d2 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'Phase', cols = ccSeurat_colors)  + ggtitle('Phase')
  d3 = DimPlot(seurat2, reduction = 'umap', label=F, group.by = 'ccAFv2', cols = ccAFv2_colors[sub2]) + ggtitle('ccAFv2')
  pdf(file.path(resdir2, paste0(tag, '.pdf')), width = 10, height = 8)
  lst = list(d1, d2, d3)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
  # gene sets
  p1 = lapply(ensembl_features1_plot, function(goi) {
    if(goi %in% rownames(seurat2)){
      fp1 = FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE, label=F, pt.size = 0.25) + ggtitle(rownames(ensembl_features1)[ensembl_features1$ensembl_features1 == goi]) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.text=element_text(size=8))
      return(fp1 + FontSize(x.title = 10, y.title = 10))
    }
  })
  grid.arrange(grobs = p1, layout_matrix = rbind(c(1, 2, 3, 4), c(5,6,7,8), c(9,10,11, 12)), top = '')
  #------- ccAFv2 vs. seurat clusters - cell percentages -------#
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(table(seurat2$ccAFv2, seurat2$seurat_clusters), beside = FALSE, col = ccAFv2_colors[sub2], xlab = 'clusters ID', ylab = 'Cell count', legend.text = rownames(table(seurat2$ccAFv2, seurat2$seurat_clusters)), args.legend=list(title='ccAFv2 classification'))
  #--- ccAFv2 vs. cluster ids stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$seurat_clusters)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$seurat_clusters)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$seurat_clusters)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 4.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 vs. ccSeurat stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$Phase)
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$Phase)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$Phase)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 1.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 and number of marker genes boxplot ---#
  v1 = VlnPlot(seurat2, features = 'ccAFv2_mgene_counts', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('ccAFv2 marker gene counts')
  v2 = VlnPlot(seurat2, features = 'nCount_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nCount_RNA')
  v3 = VlnPlot(seurat2, features = 'nFeature_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub4]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nFeature_RNA')
  lst2 = list(v2, v3, v1)
  grid.arrange(grobs = lst2, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
  print(DoHeatmap(object = seurat2, features = names(top10_symbol), group.colors = ccAFv2_colors[sub4], size = 4) + scale_y_discrete(labels = top10_symbol))
  dev.off()
  # Change factored metadata to characters
  seurat2$ccAFv2 = as.character(seurat2$ccAFv2)
  seurat2$Phase = as.character(seurat2$Phase)
  cat('saving processed data as loom and rds...\n')
  #data_loom_2 <- as.loom(seurat2, file.path(resdir3, paste0(tag, '_processed.loom')), verbose = FALSE, overwrite = TRUE)
  #data_loom_2$close_all()
  saveRDS(seurat2, file.path(resdir3, paste0(tag, '_processed.rds')))
  return(seurat2)
}

#---------------------------------------------------
# Downstream analysis & plotting (all together)
#--------------------------------------------------

scProcessData(res_dir = 'data/GSC', tag = 'BT322', resolution = 0.6)
scProcessData(res_dir = 'data/GSC', tag = 'BT324', resolution = 0.6, clusters_to_remove = c(2, 8, 10, 11, 13), norm_regress = T)
