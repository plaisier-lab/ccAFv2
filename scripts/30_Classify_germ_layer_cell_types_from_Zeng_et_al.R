##########################################################
## ccAFv2:  Classify Zeng et al. cell types             ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/testData/GSE155121:/files' cplaisier/ccafv2_seurat5

devtools::install_github("immunogenomics/presto")
library(presto)
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
library(ccAFv2)
library(scales)
library(pheatmap)
use_python('/usr/bin/python3')
remotes::install_version("matrixStats", version="1.1.0")
library(matrixStats)

scQC = function(data, tag, week_stage, v1 = 3000, v2 = 40000, h1 = 0.01, h2 = 0.12){
  seurat1 = data
  tag1 = tag
  ws1 = week_stage
  cat('\n',week_stage,'\n')
  # Create new folders to save results
  resdir3 = file.path(resdir, paste0(tag, '/', ws1))
  resdir4 = file.path(resdir3, 'analysis_output')
  resdir5 = file.path(resdir3, 'seurat_objects')
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
}

scQC_and_norm = function(data, tag, week_stage, v1 = 3000, v2 = 40000, h1 = 0.01, h2 = 0.12){
  seurat1 = data
  tag1 = tag
  ws1 = week_stage
  cat('\n',week_stage,'\n')
  # Create new folders to save results
  resdir3 = file.path(resdir, paste0(tag, '/', ws1))
  resdir4 = file.path(resdir3, 'analysis_output')
  resdir5 = file.path(resdir3, 'seurat_objects')
  seurat_subset = subset(seurat1, subset = week_stage == ws1)
  cat('\n', dim(seurat_subset),'\n')
  # Quality control filtering
  keep.detect = which(seurat_subset@meta.data$percent.mito < h2 & seurat_subset@meta.data$percent.mito > h1 & seurat_subset@meta.data$nCount_RNA < v2 & seurat_subset@meta.data$nCount_RNA > v1)
  seurat_subset = subset(seurat_subset, cells=colnames(seurat_subset)[keep.detect])
  cat('\n', dim(seurat_subset),'\n')
  cat('Normalization \n')
  seurat_subset2 = SCTransform(seurat_subset, verbose = FALSE)
  seurat_subset2 = CellCycleScoring(object=seurat_subset2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
  phase_calls = seurat_subset2$Phase
  write.csv(phase_calls, file.path(resdir4, paste0(tag, '_ccSeurat_calls.csv')))
  # Change to Ensembl ID
  seurat_subset@assays$RNA@data@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
  seurat_subset@assays$RNA@counts@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
  seurat_subset$Phase = phase_calls
  seurat_subset = SCTransform(seurat_subset, verbose = FALSE, return.only.var.genes = FALSE)
  seurat_subset = PredictCellCycle(seurat_subset, do_sctransform=FALSE)
  write.csv(seurat_subset$ccAFv2, file.path(resdir4, paste0(tag, '_ccAFv2_calls.csv')))
  tmp = data.frame(table(seurat_subset$ccAFv2))
  rownames(tmp) = tmp$Var1
  write.csv(data.frame((tmp['Freq']/dim(seurat_subset)[2])*100), file.path(resdir4, paste0(tag, '_ccAFv2_call_frequency.csv')))
  seurat_subset$ccAFv2 = as.character(seurat_subset$ccAFv2)
  seurat_subset$Main_cell_type = as.character(seurat_subset$Main_cell_type)
  seurat_subset$week_stage = as.character(seurat_subset$week_stage)
  saveRDS(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.rds')))
  #SaveLoom(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.loom')), verbose = FALSE, overwrite = TRUE)
  #SaveH5Seurat(seurat_subset, file.path(resdir5, paste0(ws1, '_normalized_ensembl.h5Seurat')), overwrite = TRUE)
  #Convert(file.path(resdir5, paste0(ws1, '_normalized_ensembl.h5Seurat')), dest = 'h5ad', overwrite = TRUE)
}

plotHeatmap = function(data, tag, save_dir = file.path(resdir2)){
  datas = data
  savedir = save_dir
  df = table(datas$seurat_clusters, datas$ccAFv2)
  enrich = data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
  rownames(enrich) = rownames(df)
  colnames(enrich) = colnames(df)
  N = sum(df)
  for (class1 in colnames(df)){
    k = sum(df[,class1])
    for (clust1 in rownames(df)){
      m = sum(df[clust1,])
      q = df[clust1, class1]
      n = N-m
      enrich[clust1,class1] = phyper(q, m, n, k, lower.tail = F)
    }
  }
  pm = -log10(as.matrix(enrich))
  pm[pm>20] = 20
  pdf(file.path(savedir, paste0(tag, '_ccAFv2_seurat_clusters_heatmap.pdf')))
  print(pheatmap(pm, cluster_cols = F, cluster_rows = F, colorRampPalette(c("white", "red"))(100), display_numbers = round(-log10(as.matrix(enrich)),2)))
  dev.off()
}

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
  seurat2 = readRDS(file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id),'.rds')))
  #nFeature_RNA <- colSums(seurat1[["RNA"]]$counts >0)  # Order ccSeurat calls
  #seurat2$nFeature_RNA = nFeature_RNA
  sub1 = ccSeurat_order %in% factor(seurat2$Phase)
  seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub1])
  # Order ccAFv2 calls
  sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
  seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
  seurat2 = RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 = RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  # Plot heatmap of hypergeometrics (cells in a ccAFv2 state cells vs. cells in a cluster)
  res1 = resolution
  plotHeatmap(seurat2, tag = res1, resdir2) # saved out as a pdf
  # Find cluster marker genes for each ccAFv2 class
  Idents(seurat2) = seurat2$ccAFv2
  cluster_markers= FindAllMarkers(seurat2, only.pos = TRUE, logfc.threshold = 0.25)
  cluster_markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 0.5) %>%
      slice_head(n = 10) %>%
      ungroup() -> top10
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
  #--- ccAFv2 vs. cluster ids stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$seurat_clusters)
  cf <- cf[rowSums(cf[])>0,]
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=ncol(cf), nrow=nrow(cf))
  for(i in c(1:ncol(cf))){
    for(n in c(1:nrow(cf))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[nrow(cf)+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 4.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 vs. ccSeurat stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$Phase)
  cf <- cf[rowSums(cf[])>0,]
  totals <- colSums(cf)
  data.frame(totals)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=ncol(cf), nrow=nrow(cf))
  for(i in c(1:ncol(cf))){
    for(n in c(1:nrow(cf))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[nrow(cf)+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cf)
  sub4 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(8, 8, 8, 8) + 2.0)
  barplot(cf_1, xlab = '', ylab = 'Cell Percentage', las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub4], args.legend=list(x=ncol(cf_1) + 1.5, y=max(colSums(cf_1)), bty = 'n'))
  #--- ccAFv2 and number of marker genes boxplot ---#
  #v2 = VlnPlot(seurat2, features = 'nCount_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub2]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nCount_RNA')
  #v3 = VlnPlot(seurat2, features = 'nFeature_RNA', group.by = 'ccAFv2', cols = ccAFv2_colors[sub2]) + theme(legend.position = 'none') + xlab('ccAFv2') + ylab('nFeature_RNA')
  #lst2 = list(v2, v3)
  #grid.arrange(grobs = lst2, layout_matrix = rbind(c(1, 2), c(3, NA)), top = '')
  print(DoHeatmap(object = seurat2, features = names(top10_symbol), size = 4, group.colors = ccAFv2_colors[sub2]) + scale_y_discrete(labels = top10_symbol))
  print(RidgePlot(seurat2, features = ensembl_features3_plot, ncol=2, cols = ccAFv2_colors[sub2]))
  dev.off()
  # Change factored metadata to characters
  seurat2$ccAFv2 = as.character(seurat2$ccAFv2)
  seurat2$Phase = as.character(seurat2$Phase)
  #cat('saving processed data as loom and rds...\n')
  #data_loom_2 <- as.loom(seurat2, file.path(resdir3, paste0(tag, '_processed.loom')), verbose = FALSE, overwrite = TRUE)
  #data_loom_2$close_all()
  saveRDS(seurat2, file.path(resdir3, paste0(tag, '_processed.rds')))
  return(seurat2)
}

plotting = list()
plotStacked = function(res_dir, tag, save_dir = 'analysis_output'){
  cat('\n',tag,'\n')
  # Set up folders
  resdir1 = file.path(resdir, res_dir, tag)
  resdir2 = file.path(resdir1, save_dir)
  #---------------------------
  # Load ccAFv2 csv
  #---------------------------
  freq = read.csv(file.path(resdir2, paste0(res_dir, '_ccAFv2_call_frequency.csv')), row.names = 'X')
  plotting[[paste0(res_dir, '_', tag)]] = freq
  sub2 = ccAFv2_order %in% rownames(freq)
  freq$ccAFv2 <- rownames(freq)
  pdf(file.path(resdir, paste0(res_dir, '_', tag, '_ccAFv2_stacked_bar_plot.pdf')))
  ggplot(freq, aes(x = "", y = Freq, fill = ccAFv2)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(x = NULL, y = NULL, fill = "Category") +
    theme_void() +
    theme(legend.position = "right")
  dev.off()
}


all_mgenes = list()
find_G0genes = function(res_dir, tag){
  cat('\n',tag,'\n')
  # Set up folders
  resdir1 = file.path(resdir, res_dir, tag)
  resdir2 = file.path(resdir1, 'analysis_output')
  mgenes = read.csv(file.path(resdir2, paste0(tag, '_scTransform_Markers_together.csv')), row.names = 'X')
  mgenes2 = mgenes[mgenes$cluster == 'Neural G0',]
  mgenes2[mgenes2['p_val_adj'] <= 0.05,]
  all_mgenes[[paste0(resdir, '_', tag)]] = rownames(mgenes2)
  return(all_mgenes)
}

find_G0genes(res_dir = 'Endoderm', tag = 'W4-1')

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

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

features3 <- c("CCND1", "CCNE2", "CCNA2", "CCNB1", "CDK1", "CDK2")
# convert to ensembl IDs
ensembl_features3 = mapIds(org.Hs.eg.db, keys = features3, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features3 = na.omit(data.frame(ensembl_features3))
ensembl_features3_plot = ensembl_features3$ensembl_features3

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccAFv2_order_opp = c('Unknown', 'M/Early G1', 'G2/M', 'S/G2', 'S', 'Late G1', 'G1', 'Neural G0')

# Set up directory structure
#setwd("files/")
resdir = 'new'

for(celltype1 in c('Notochord', 'Endoderm', 'IPC', 'Endothelium', 'Blood', 'NMP', 'Erythroid', 'Immune Cell', 'Cholin N', 'Epithelium', 'LPM', 'Neural Crest', 'Osteoblast', 'Chondroblast', 'Glu N', 'Myoblast', 'DA N', 'GABA N', 'MesoPro')){
  resdir2 = file.path(resdir)
  #Convert(file.path(resdir2, paste0(celltype1, '.h5ad')), dest = 'h5seurat', overwrite = TRUE)
  seurat1 <- LoadH5Seurat(file.path(resdir2,paste0(celltype1,'.h5seurat')))
  mito.genes <- grep("MT-", rownames(seurat1))
  percent.mito = colSums(seurat1[mito.genes,])/colSums(seurat1@assays[['RNA']])
  seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")
  nCount = colSums(x = seurat1, slot = "counts")  # nCount_RNA
  seurat1 <- AddMetaData(seurat1, nCount, col.name = "nCount_RNA")
  # Create folder for each main cell type
  dir.create(file.path(resdir2, celltype1), showWarnings = FALSE)
  # Create new folders to save results
  resdir3 = file.path(resdir2, celltype1)
  saveRDS(seurat1, file.path(resdir3, paste0(celltype1, '_unfiltered.rds')))
}

celltype1 = 'Notochord'
resdir3 = file.path(resdir, celltype1)
seurat1 = readRDS(file.path(resdir3, paste0(celltype1, '_unfiltered.rds')))
table(seurat1$week_stage)
table(seurat1$cell_type)


# Notochord
scQC(seurat1, 'Notochord', 'W5-2', v2 = 16000)
scQC_and_norm(seurat1, 'Notochord', 'W5-2', v2 = 16000)
scProcessData(res_dir = 'new/Notochord', tag = 'W5-2')

# Endoderm
scQC(seurat1, 'Endoderm', 'W4-1', h2 = 0.10, v2 = 40000)
scQC_and_norm(seurat1, 'Endoderm', 'W4-1', h2 = 0.10, v2 = 40000)
scProcessData(res_dir = 'new/Endoderm', tag = 'W4-1', resolution = 0.3)
plotStacked(res_dir = 'Endoderm', tag = 'W4-1')

# IPC
scQC(seurat1, 'IPC', 'W9-1', h1 = 0.005, h2 = 0.025, v2 = 22000)
scQC(seurat1, 'IPC', 'W9-2', h1 = 0.005, h2 = 0.05, v1 = 4000, v2 = 24000)
scQC(seurat1, 'IPC', 'W12-1', h1 = 0.005, h2 = 0.03, v1 = 4000, v2 = 19000)
scQC_and_norm(seurat1, 'IPC', 'W9-1', h1 = 0.005, h2 = 0.025, v2 = 22000)
scProcessData(res_dir = 'new/IPC', tag = 'W9-1', resolution = 0.6)
scQC_and_norm(seurat1, 'IPC', 'W9-2', h1 = 0.005, h2 = 0.05, v1 = 4000, v2 = 24000)
scProcessData(res_dir = 'new/IPC', tag = 'W9-2')
scQC_and_norm(seurat1, 'IPC', 'W12-1', h1 = 0.005, h2 = 0.03, v1 = 4000, v2 = 19000)
scProcessData(res_dir = 'new/IPC', tag = 'W12-1', resolution = 0.7)

# Endothelium
scQC(seurat1, 'Endothelium', 'W4-1')
scQC(seurat1, 'Endothelium', 'W4-3', h2 = 0.06, v2 = 35000)
scQC_and_norm(seurat1, 'Endothelium', 'W4-1')
scQC_and_norm(seurat1, 'Endothelium', 'W4-3', h2 = 0.06, v2 = 35000)
scProcessData(res_dir = 'new/Endothelium', tag = 'W4-1')
scProcessData(res_dir = 'new/Endothelium', tag = 'W4-3')

# Blood
scQC(seurat1, 'Blood', 'W9-2', h1 = 0, h2 = 0.03, v1 = 7000)
scQC(seurat1, 'Blood', 'W12-1', h1 = 0, h2 = 0.03, v1 = 7000, v2 = 30000)
scQC_and_norm(seurat1, 'Blood', 'W9-2', h1 = 0, h2 = 0.03, v1 = 7000)
scQC_and_norm(seurat1, 'Blood', 'W12-1', h1 = 0, h2 = 0.03, v1 = 7000, v2 = 30000)
scProcessData(res_dir = 'new/Blood', tag = 'W9-2')
scProcessData(res_dir = 'new/Blood', tag = 'W12-1')

# NMP
scQC(seurat1, 'NMP', 'W3-1', h1 = 0, h2 = 0.025, v1 = 8000)
scQC_and_norm(seurat1, 'NMP', 'W3-1', h1 = 0, h2 = 0.025, v1 = 8000)
scProcessData(res_dir = 'new/NMP', tag = 'W3-1', resolution = 0.6)

# Erythroid
scQC(seurat1, 'Erythroid', 'W4-2', h1 = 0, h2 = 0.025, v1 = 8000)
scQC(seurat1, 'Erythroid', 'W6-1', h1 = 0, h2 = 0.015, v1 = 8000, v2 = 60000)
scQC_and_norm(seurat1, 'Erythroid', 'W4-2', h1 = 0, h2 = 0.025, v1 = 8000)
scQC_and_norm(seurat1, 'Erythroid', 'W6-1', h1 = 0, h2 = 0.015, v1 = 8000, v2 = 60000)
scProcessData(res_dir = 'new/Erythroid', tag = 'W4-2')
scProcessData(res_dir = 'new/Erythroid', tag = 'W6-1')

# Immune cells
scQC(seurat1, 'Immune Cell', 'W9-2', h2 = 0.10)
scQC(seurat1, 'Immune Cell', 'W12-1', h2 = 0.06, v2 = 20000)
scQC_and_norm(seurat1, 'Immune Cell', 'W9-2', h2 = 0.10)
scQC_and_norm(seurat1, 'Immune Cell', 'W12-1', h2 = 0.06, v2 = 20000)
scProcessData(res_dir = 'new/Immune Cell', tag = 'W9-2')
scProcessData(res_dir = 'new/Immune Cell', tag = 'W12-1')

# Cholin N
scQC(seurat1, 'Cholin N', 'W5-1', h2 = 0.04, v1 = 6000, v2 = 23000)
scQC_and_norm(seurat1, 'Cholin N', 'W5-1', h2 = 0.04, v1 = 6000, v2 = 23000)
scProcessData(res_dir = 'new/Cholin N', 'W5-1')

# Epithelium
scQC(seurat1, 'Epithelium', 'W4-1', v1 = 5000)
scQC(seurat1, 'Epithelium', 'W5-2', h1 = 0.03, v2 = 26000)
scQC_and_norm(seurat1, 'Epithelium', 'W4-1', v1 = 5000)
scQC_and_norm(seurat1, 'Epithelium', 'W5-2', h1 = 0.03, v2 = 26000)
scProcessData(res_dir = 'new/Epithelium', 'W4-1', resolution = 0.6)
scProcessData(res_dir = 'new/Epithelium', 'W5-2', resolution = 0.6)

# LPM
scQC(seurat1, 'LPM', 'W3-1', h1 = 0, h2 = 0.03, v1 = 7000)
scQC(seurat1, 'LPM', 'W4-1', v1 = 5000, v2 = 32000)
scQC(seurat1, 'LPM', 'W5-1', v1 = 5000, v2 = 25000, h2 = 0.06)
scQC_and_norm(seurat1, 'LPM', 'W3-1', h1 = 0, h2 = 0.03, v1 = 7000)
scQC_and_norm(seurat1, 'LPM', 'W4-1', v1 = 5000, v2 = 32000)
scQC_and_norm(seurat1, 'LPM', 'W5-1', v1 = 5000, v2 = 25000, h2 = 0.06)
scProcessData(res_dir = 'new/LPM', 'W3-1', resolution = 0.6)
scProcessData(res_dir = 'new/LPM', 'W4-1', resolution = 0.5)
scProcessData(res_dir = 'new/LPM', 'W5-1')

# Neural Crest
scQC(seurat1, 'Neural Crest', 'W4-1', v1 = 5000, v2 = 30000, h2 = 0.10)
scQC(seurat1, 'Neural Crest', 'W4-2', v1 = 5000, v2 = 30000, h2 = 0.05)
scQC(seurat1, 'Neural Crest', 'W4-3', h2 = 0.04, v2 = 17000)
scQC(seurat1, 'Neural Crest', 'W5-1', v1 = 5000, h2 = 0.05)
scQC(seurat1, 'Neural Crest', 'W5-2', v2 = 26000)
scQC(seurat1, 'Neural Crest', 'W5-3', v1 = 6000, h2 = 0.06)
scQC_and_norm(seurat1, 'Neural Crest', 'W4-1', v1 = 5000, v2 = 30000, h2 = 0.10)
scQC_and_norm(seurat1, 'Neural Crest', 'W4-2', v1 = 5000, v2 = 30000, h2 = 0.05)
scQC_and_norm(seurat1, 'Neural Crest', 'W4-3', h2 = 0.04, v2 = 17000)
scQC_and_norm(seurat1, 'Neural Crest', 'W5-1', v1 = 5000, h2 = 0.05)
scQC_and_norm(seurat1, 'Neural Crest', 'W5-2', v2 = 26000)
scQC_and_norm(seurat1, 'Neural Crest', 'W5-3', v1 = 6000, h2 = 0.06)
scProcessData(res_dir = 'new/Neural Crest', 'W4-1')
scProcessData(res_dir = 'new/Neural Crest', 'W4-2')
scProcessData(res_dir = 'new/Neural Crest', 'W4-3')
scProcessData(res_dir = 'new/Neural Crest', 'W5-1', resolution = 0.6)
scProcessData(res_dir = 'new/Neural Crest', 'W5-2')
scProcessData(res_dir = 'new/Neural Crest', 'W5-3')

# Osteoblast
scQC(seurat1, 'Osteoblast', 'W4-2', h2 = 0.05, v2 = 15000)
scQC(seurat1, 'Osteoblast', 'W5-2', v2 = 15000)
scQC(seurat1, 'Osteoblast', 'W6-1', v1 = 5000, v2 = 20000, h2 = 0.04)
scQC(seurat1, 'Osteoblast', 'W7-1', h2 = 0.04, v2 = 15000)
scQC(seurat1, 'Osteoblast', 'W8-1', h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'Osteoblast', 'W4-2', h2 = 0.05, v2 = 15000)
scQC_and_norm(seurat1, 'Osteoblast', 'W5-2', v2 = 15000)
scQC_and_norm(seurat1, 'Osteoblast', 'W6-1', v1 = 5000, v2 = 20000, h2 = 0.04)
scQC_and_norm(seurat1, 'Osteoblast', 'W7-1', h2 = 0.04, v2 = 15000)
scQC_and_norm(seurat1, 'Osteoblast', 'W8-1', h2 = 0.04, v2 = 20000)
scProcessData(res_dir = 'new/Osteoblast', 'W4-2')
scProcessData(res_dir = 'new/Osteoblast', 'W5-2', resolution = 0.6)
scProcessData(res_dir = 'new/Osteoblast', 'W6-1')
scProcessData(res_dir = 'new/Osteoblast', 'W7-1')
scProcessData(res_dir = 'new/Osteoblast', 'W8-1', resolution = 0.6)

# Chondroblast
scQC(seurat1, 'Chondroblast', 'W4-2', h2 = 0.05, v2 = 20000)
scQC(seurat1, 'Chondroblast', 'W4-3', h2 = 0.05, v2 = 15000)
scQC(seurat1, 'Chondroblast', 'W5-1', v1 = 5000, h2 = 0.05, v2 = 23000)
scQC(seurat1, 'Chondroblast', 'W5-2', v2 = 20000)
scQC(seurat1, 'Chondroblast', 'W5-3', v1 = 5000, v2 = 27000, h2 = 0.05)
scQC(seurat1, 'Chondroblast', 'W6-1', h2 = 0.035, v1 = 5000, v2 = 27000)
scQC(seurat1, 'Chondroblast', 'W7-1', h2 = 0.04, v2 = 20000)
scQC(seurat1, 'Chondroblast', 'W8-1', h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'Chondroblast', 'W4-2', h2 = 0.05, v2 = 20000)
scQC_and_norm(seurat1, 'Chondroblast', 'W4-3', h2 = 0.05, v2 = 15000)
scQC_and_norm(seurat1, 'Chondroblast', 'W5-1', v1 = 5000, h2 = 0.05, v2 = 23000)
scQC_and_norm(seurat1, 'Chondroblast', 'W5-2', v2 = 20000)
scQC_and_norm(seurat1, 'Chondroblast', 'W5-3', v1 = 5000, v2 = 27000, h2 = 0.05)
scQC_and_norm(seurat1, 'Chondroblast', 'W6-1', h2 = 0.035, v1 = 5000, v2 = 27000)
scQC_and_norm(seurat1, 'Chondroblast', 'W7-1', h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'Chondroblast', 'W8-1', h2 = 0.04, v2 = 20000)
scProcessData(res_dir = 'new/Chondroblast', 'W4-2')
scProcessData(res_dir = 'new/Chondroblast', 'W4-3')
scProcessData(res_dir = 'new/Chondroblast', 'W5-1', resolution = 0.6)
scProcessData(res_dir = 'new/Chondroblast', 'W5-2', resolution = 0.6)
scProcessData(res_dir = 'new/Chondroblast', 'W5-3')
scProcessData(res_dir = 'new/Chondroblast', 'W6-1')
scProcessData(res_dir = 'new/Chondroblast', 'W7-1')
scProcessData(res_dir = 'new/Chondroblast', 'W8-1')

# Glu N
scQC(seurat1, 'Glu N', 'W4-2', v1 = 5000, v2 = 27000, h2 = 0.04)
scQC(seurat1, 'Glu N', 'W5-1', v1 = 10000, h2 = 0.04)
scQC(seurat1, 'Glu N', 'W5-2', v1 = 5000, v2 = 25000, h2 = 0.10)
scQC(seurat1, 'Glu N', 'W6-1', h1 = 0, h2 = 0.04, v1 = 5000, v2 = 35000)
scQC(seurat1, 'Glu N', 'W9-1', h1 = 0.005, h2 = 0.03, v2 = 16000)
scQC(seurat1, 'Glu N', 'W9-2', h1 = 0.005, h2 = 0.04, v1 = 5000, v2 = 23000)
scQC(seurat1, 'Glu N', 'W12-1', h1 = 0.005, h2 = 0.035, v2 = 17000)
scQC_and_norm(seurat1, 'Glu N', 'W4-2', v1 = 5000, v2 = 27000, h2 = 0.04)
scQC_and_norm(seurat1, 'Glu N', 'W5-1', v1 = 10000, h2 = 0.04)
scQC_and_norm(seurat1, 'Glu N', 'W5-2', v1 = 5000, v2 = 25000, h2 = 0.10)
scQC_and_norm(seurat1, 'Glu N', 'W6-1', h1 = 0, h2 = 0.04, v1 = 5000, v2 = 35000)
scQC_and_norm(seurat1, 'Glu N', 'W9-1', h1 = 0.005, h2 = 0.03, v2 = 16000)
scQC_and_norm(seurat1, 'Glu N', 'W12-1', h1 = 0.005, h2 = 0.035, v2 = 17000)
scProcessData(res_dir = 'new/Glu N', 'W4-2')
scProcessData(res_dir = 'new/Glu N', 'W5-1')
scProcessData(res_dir = 'new/Glu N', 'W5-2')
scProcessData(res_dir = 'new/Glu N', 'W6-1')
scProcessData(res_dir = 'new/Glu N', 'W9-1')
scProcessData(res_dir = 'new/Glu N', 'W12-1')

# Myoblast
scQC(seurat1, 'Myoblast', 'W3-1', h1 = 0, h2 = 0.03, v1 = 5000, v2 = 30000)
scQC(seurat1, 'Myoblast', 'W4-1', v1 = 5000, v2 = 30000)
scQC(seurat1, 'Myoblast', 'W4-2', h2 = 0.06, v2 = 20000)
scQC(seurat1, 'Myoblast', 'W4-3', v2 = 17000, h2 = 0.05)
scQC(seurat1, 'Myoblast', 'W5-1', h2 = 0.06, v1 = 5000, v2 = 30000)
scQC(seurat1, 'Myoblast', 'W5-2', v2 = 25000)
scQC(seurat1, 'Myoblast', 'W5-3', v1 = 5000, v2 = 30000, h2 = 0.05)
scQC(seurat1, 'Myoblast', 'W6-1', h2 = 0.04, v2 = 23000, v1 = 5000)
scQC(seurat1, 'Myoblast', 'W7-1', h2 = 0.04, v2 = 17000)
scQC_and_norm(seurat1, 'Myoblast', 'W3-1', h1 = 0, h2 = 0.03, v1 = 5000, v2 = 30000)
scQC_and_norm(seurat1, 'Myoblast', 'W4-1', v1 = 5000, v2 = 30000)
scQC_and_norm(seurat1, 'Myoblast', 'W4-2', h2 = 0.06, v2 = 20000)
scQC_and_norm(seurat1, 'Myoblast', 'W4-3', v2 = 17000, h2 = 0.05)
scQC_and_norm(seurat1, 'Myoblast', 'W5-1', h2 = 0.06, v1 = 5000, v2 = 30000)
scQC_and_norm(seurat1, 'Myoblast', 'W5-2', v2 = 25000)
scQC_and_norm(seurat1, 'Myoblast', 'W5-3', v1 = 5000, v2 = 30000, h2 = 0.05)
scQC_and_norm(seurat1, 'Myoblast', 'W6-1', h2 = 0.04, v2 = 23000, v1 = 5000)
scQC_and_norm(seurat1, 'Myoblast', 'W7-1', h2 = 0.04, v2 = 17000)
scProcessData(res_dir = 'new/Myoblast', tag = 'W3-1')
scProcessData(res_dir = 'new/Myoblast', tag = 'W4-1')
scProcessData(res_dir = 'new/Myoblast', tag = 'W4-2')
scProcessData(res_dir = 'new/Myoblast', tag = 'W4-3')
scProcessData(res_dir = 'new/Myoblast', tag = 'W5-1')
scProcessData(res_dir = 'new/Myoblast', tag = 'W5-2')
scProcessData(res_dir = 'new/Myoblast', tag = 'W5-3')
scProcessData(res_dir = 'new/Myoblast', tag = 'W6-1')
scProcessData(res_dir = 'new/Myoblast', tag = 'W7-1')

# DA N
scQC(seurat1, 'DA N', 'W4-2', h1 = 0, h2 = 0.04, v2 = 23000)
scQC(seurat1, 'DA N', 'W5-1', v1 = 5000, v2 = 25000, h2 = 0.04, h1 = 0.005)
scQC(seurat1, 'DA N', 'W5-2', v2 = 20000, h2 = 0.1)
scQC(seurat1, 'DA N', 'W6-1', h1 = 0.005, h2 = 0.03, v1 = 5000, v2 = 28000)
scQC(seurat1, 'DA N', 'W7-1', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC(seurat1, 'DA N', 'W8-1', h2 = 0.035, v2 = 25000, h1 = 0.005)
scQC(seurat1, 'DA N', 'W9-1', h1 = 0.005, h2 = 0.05, v2 = 15000)
scQC(seurat1, 'DA N', 'W9-2', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC(seurat1, 'DA N', 'W12-1', h1 = 0.005, h2 = 0.04, v2 = 15000)
scQC_and_norm(seurat1, 'DA N', 'W4-2', h1 = 0, h2 = 0.04, v2 = 23000)
scQC_and_norm(seurat1, 'DA N', 'W5-1', v1 = 5000, v2 = 25000, h2 = 0.04, h1 = 0.005)
scQC_and_norm(seurat1, 'DA N', 'W5-2', v2 = 20000, h2 = 0.1)
scQC_and_norm(seurat1, 'DA N', 'W6-1', h1 = 0.005, h2 = 0.03, v1 = 5000, v2 = 28000)
scQC_and_norm(seurat1, 'DA N', 'W7-1', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'DA N', 'W8-1', h2 = 0.035, v2 = 25000, h1 = 0.005)
scQC_and_norm(seurat1, 'DA N', 'W9-1', h1 = 0.005, h2 = 0.05, v2 = 15000)
scQC_and_norm(seurat1, 'DA N', 'W9-2', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'DA N', 'W12-1', h1 = 0.005, h2 = 0.04, v2 = 15000)
scProcessData(res_dir = 'new/DA N', tag = 'W4-2')
scProcessData(res_dir = 'new/DA N', tag = 'W5-1')
scProcessData(res_dir = 'new/DA N', tag = 'W5-2')
scProcessData(res_dir = 'new/DA N', tag = 'W6-1', resolution = 0.6)
scProcessData(res_dir = 'new/DA N', tag = 'W7-1')
scProcessData(res_dir = 'new/DA N', tag = 'W8-1')
scProcessData(res_dir = 'new/DA N', tag = 'W9-1')
scProcessData(res_dir = 'new/DA N', tag = 'W9-2')
scProcessData(res_dir = 'new/DA N', tag = 'W12-1')

# GABA N
scQC(seurat1, 'GABA N', 'W4-2', h2 = 0.05, v2 = 25000)
scQC(seurat1, 'GABA N', 'W6-1', h1 = 0.005, h2 = 0.03, v1 = 5000, v2 = 30000)
scQC(seurat1, 'GABA N', 'W7-1', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC(seurat1, 'GABA N', 'W8-1', h1 = 0.005, h2 = 0.04, v2 = 23000)
scQC(seurat1, 'GABA N', 'W9-1', h1 = 0.005, h2 = 0.04, v2 = 18000)
scQC(seurat1, 'GABA N', 'W9-2', h2 = 0.05, v2 = 25000)
scQC(seurat1, 'GABA N', 'W12-1', h1 = 0.005, h2 = 0.04, v2 = 17000)
scQC_and_norm(seurat1, 'GABA N', 'W4-2', h2 = 0.05, v2 = 25000)
scQC_and_norm(seurat1, 'GABA N', 'W6-1', h1 = 0.005, h2 = 0.03, v1 = 5000, v2 = 30000)
scQC_and_norm(seurat1, 'GABA N', 'W7-1', h1 = 0.005, h2 = 0.04, v2 = 20000)
scQC_and_norm(seurat1, 'GABA N', 'W8-1', h1 = 0.005, h2 = 0.04, v2 = 23000)
scQC_and_norm(seurat1, 'GABA N', 'W9-1', h1 = 0.005, h2 = 0.04, v2 = 18000)
scQC_and_norm(seurat1, 'GABA N', 'W9-2', h2 = 0.05, v2 = 25000)
scQC_and_norm(seurat1, 'GABA N', 'W12-1', h1 = 0.005, h2 = 0.04, v2 = 17000)
scProcessData(res_dir = 'new/GABA N', tag = 'W4-2')
scProcessData(res_dir = 'new/GABA N', tag = 'W6-1')
scProcessData(res_dir = 'new/GABA N', tag = 'W7-1')
scProcessData(res_dir = 'new/GABA N', tag = 'W8-1')
scProcessData(res_dir = 'new/GABA N', tag = 'W9-1')
scProcessData(res_dir = 'new/GABA N', tag = 'W9-2')
scProcessData(res_dir = 'new/GABA N', tag = 'W12-1')

# MesoPro
scQC(seurat1, 'MesoPro', 'W3-1', h1 = 0, h2 = 0.03, v1 = 5000)
scQC(seurat1, 'MesoPro', 'W4-1', v1 = 5000, v2 = 38000)
scQC(seurat1, 'MesoPro', 'W4-2', h2 = 0.06, v2 = 30000)
scQC(seurat1, 'MesoPro', 'W4-3', h2 = 0.06, v2 = 23000)
scQC(seurat1, 'MesoPro', 'W5-1', v1 = 5000, h2 = 0.05)
scQC(seurat1, 'MesoPro', 'W5-2', v1 = 5000, v2 = 30000, h1 = 0.02)
scQC(seurat1, 'MesoPro', 'W5-3', v1 = 5000, h2 = 0.05)
scQC(seurat1, 'MesoPro', 'W6-1', h2 = 0.04, v1 = 5000, v2 = 38000)
scQC(seurat1, 'MesoPro', 'W7-1', h2 = 0.04, v1 = 5000, v2 = 25000)
scQC(seurat1, 'MesoPro', 'W8-1', h2 = 0.05, v2 = 25000)
scQC_and_norm(seurat1, 'MesoPro', 'W3-1', h1 = 0, h2 = 0.03, v1 = 5000)
scQC_and_norm(seurat1, 'MesoPro', 'W4-1', v1 = 5000, v2 = 38000)
scQC_and_norm(seurat1, 'MesoPro', 'W4-2', h2 = 0.06, v2 = 30000)
scQC_and_norm(seurat1, 'MesoPro', 'W4-3', h2 = 0.06, v2 = 23000)
scQC_and_norm(seurat1, 'MesoPro', 'W5-1', v1 = 5000, h2 = 0.05)
scQC_and_norm(seurat1, 'MesoPro', 'W5-2', v1 = 5000, v2 = 30000, h1 = 0.02)
scQC_and_norm(seurat1, 'MesoPro', 'W5-3', v1 = 5000, h2 = 0.05)
scQC_and_norm(seurat1, 'MesoPro', 'W6-1', h2 = 0.04, v1 = 5000, v2 = 38000)
scQC_and_norm(seurat1, 'MesoPro', 'W7-1', h2 = 0.04, v1 = 5000, v2 = 25000)
scQC_and_norm(seurat1, 'MesoPro', 'W8-1', h2 = 0.05, v2 = 25000)
scProcessData(res_dir = 'new/MesoPro', tag = 'W3-1')
scProcessData(res_dir = 'new/MesoPro', tag = 'W4-1')
scProcessData(res_dir = 'new/MesoPro', tag = 'W4-2')
scProcessData(res_dir = 'new/MesoPro', tag = 'W4-3')
scProcessData(res_dir = 'new/MesoPro', tag = 'W5-1')
scProcessData(res_dir = 'new/MesoPro', tag = 'W5-2')
scProcessData(res_dir = 'new/MesoPro', tag = 'W5-3')
scProcessData(res_dir = 'new/MesoPro', tag = 'W6-1')
scProcessData(res_dir = 'new/MesoPro', tag = 'W7-1')
scProcessData(res_dir = 'new/MesoPro', tag = 'W8-1')
