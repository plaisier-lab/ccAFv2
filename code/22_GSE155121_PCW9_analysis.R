##########################################################
## ccAFv2:  Application to GSE155121 PCW 9-1            ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2/data/GSE155121:/files' cplaisier/ccafv2_seurat4

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
use_python('/usr/bin/python3')

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Plotting order & colors
ccSeurat_order = c("G1", "S", "G2M")
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_order_opp = c('M/Early G1', 'G2/M', 'S/G2', 'S', 'Late G1', 'G1', 'Neural G0')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
tricycle_order = c('G1/G0', 'S', 'G2/M', 'M')
tricycle_colors = c("G1/G0" = "#f37f73", "S" = "#8571b2", "G2/M" = "#3db270", "M" = "#6d90ca", "Unknown" = "#D3D3D3")
ccSchwabe_order = c('G1.S', 'S', 'G2', 'G2.M', 'M.G1', 'Unknown')
ccSchwabe_colors = c("G1.S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2.M" = "#3db270", "M.G1"= "#6d90ca", "Unknown" = "#D3D3D3")
recat_order = c('G1', 'G1S', 'S', 'G2', 'G2M', 'M')
recat_colors = c("G1" = "#f37f73", "G1S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2M" = "#3db270", "M"= "#6d90ca")

ccAF_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'G1/other')
ccAF_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "G1/other" = "#9aca3c")


# cyclins to plot
features3 <- c("CCND1", "CCNE2", "CCNA2", "CCNB1", "CDK1") # original
#features3 <- c("CCND1", "CCNE2", "CCNA2", "CCNB1")
# convert to ensembl IDs
ensembl_features3 = mapIds(org.Hs.eg.db, keys = features3, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_features3 = na.omit(data.frame(ensembl_features3))
ensembl_features3_plot = ensembl_features3$ensembl_features3

neuralGO_geneset = c("SCRG1", "PLP1", "S100B", "GPM6B", "BEX1", "PTPRZ1", "PRCP", "PTN", "SOX4", "SAT1")
G1_geneset = c("IGFBP3", "IGFBP5", "MIAT", "MAP3K7CL", "AHNAK2", "TPST2", "DLG1" , "CMTM7", "C6orf15", "GJB2")
late_G1_geneset = c("EDN1", "CYR61", "ANKRD1", "CTGF", "PLK2", "UGCG", "ARID5B", "PLAU", "CCL2")
S_geneset = c("CCNE2", "CLSPN", "GINS2", "PCNA", "ATAD2", "MCM7", "MCM3", "SLBP", "GMNN", "KIAA0101")
s_g2_geneset = c("HIST1H4C", "CDK1", "HIST1H1E", "HIST1H1B", "UBE2C", "RRM2", "ZWINT", "HIST1H1C", "HMGB2")
G2_M_geneset = c("CCNB1", "CENPF", "CKS2", "PTTG1", "CDC20", "TOP2A", "NUSAP1", "CENPA")
M_early_G1_geneset = c("HMGN2", "TUBA1B", "STMN1", "BIRC5", "HMGB1", "TROAP", "HNRNPA2B1", "H2AFZ", "ARL6IP1")
goi_lst = list(neuralGO_geneset,G1_geneset, late_G1_geneset, S_geneset, s_g2_geneset, G2_M_geneset, M_early_G1_geneset )
sapply(goi_lst, length)
# convert to ensembl IDs
goi_lst_ensembl = mapIds(org.Hs.eg.db, keys = unlist(goi_lst), keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
goi_lst_ensembl2 = goi_lst_ensembl[!is.na(goi_lst_ensembl)]

mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0('human_ensembl')]

# Set up directory structure
#setwd("files/")
resdir = 'NSC'
tag1 = 'W9-1'
resdir2 = file.path(resdir, tag1)

# Read in data
seurat1 = readRDS(file.path(resdir2, 'seurat_objects/',paste0(tag1, '_normalized_ensembl.rds')))

# Order ccAFv2 cell cycle states
sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 = factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])
sub2 = ccSeurat_order %in% factor(seurat1$Phase)
seurat1$Phase <- factor(seurat1$Phase, levels = ccSeurat_order[sub2])

# Remove unknown cells
seurat_subset = subset(seurat1, subset = ccAFv2 != 'Unknown')

# Plot cyclin violin plots
pdf(file.path(resdir2, paste0('W9-1_cyclin_vln_plot_SCT_scale_data.pdf')), height = 8, width = 8)
#VlnPlot(seurat_subset, features = ensembl_features3_plot, cols = ccAFv2_colors[sub1], ncol = 2, group.by = 'ccAFv2', slot = 'data', sort = TRUE)
VlnPlot(seurat_subset, features = ensembl_features3_plot, cols = ccAFv2_colors[sub1], ncol = 2, group.by = 'ccAFv2', slot = 'scale.data')
dev.off()

# Find average expression of cyclins across cell cycle states
datas = c()
datas_std = c()
for (state1 in ccAFv2_order[sub1]){
  print(state1)
  datas[[state1]] = c()
  datas_std[[state1]] = c()
  subset1  = subset(seurat_subset, subset = ccAFv2 == state1)
  for(cyclin1 in ensembl_features3_plot){
    # Mean
    datas[[state1]][[cyclin1]] = mean(subset1[cyclin1,]$SCT@scale.data)
    # Standard deviation
    datas_std[[state1]][[cyclin1]] = sd(subset1[cyclin1,]$SCT@scale.data)
    # Median
    #datas[[state1]][[cyclin1]] = median(subset1[cyclin1,]$SCT@data)
  }
}

df_median = do.call(rbind.data.frame, datas)
colnames(df_median) = features3
df_median$ccAFv2 = rownames(df_median)
rownames(df_median) = c(1, 2, 3, 4, 5, 6, 7)

df_std = do.call(rbind.data.frame, datas_std)
colnames(df_std) = features3
df_std$ccAFv2 = rownames(df_std)
rownames(df_std) = c(1, 2, 3, 4, 5, 6, 7)

melted = melt(df_median)
melted_std = melt(df_std)
melted$std = melted_std$value
melted$ccAFv2 <- factor(melted$ccAFv2, levels=ccAFv2_order)

#pdf(file.path(resdir2, 'W9-1_cyclin_mean_expression_sct_data_ccnd2.pdf'))
#pdf(file.path(resdir2, 'W9-1_cyclin_mean_expression_sct_data_ccnd1_cdk1_082924.pdf'))
#pdf(file.path(resdir2, 'W9-1_cyclin_mean_expression_sct_data_ccnd1_with_all_cdks_082924.pdf'))
pdf(file.path(resdir2, 'W9-1_cyclin_mean_expression_sct_data_ccnd1_with_all_cdks_090324_with_std_scale_data.pdf'))
melted %>%
ggplot(aes(x = ccAFv2, y = value, group = ccAFv2, fill = ccAFv2, color = variable)) +
geom_point() +
geom_errorbar(aes(ymin=value-std, ymax=value+std), width=.2, position=position_dodge(0.05)) +
facet_wrap( ~ variable) + theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

stats = c()
for (gene1 in ensembl_features3_plot){
  print(gene1)
  stats[[gene1]] = c()
  for (state1 in c('Neural G0', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')){
     tmp = FindMarkers(seurat1, features = gene1,  group.by = 'ccAFv2', ident.2 = 'G1', ident.1 = state1, test.use = 't')
     stats[[gene1]][[state1]] = tmp$p_val_adj
  }
}

# ccAFv2 important features
important = c('PRCP', 'GPM6B', 'SCRG1','HMGN2', 'SMS', 'HMGB1', 'CCN1', 'CCN2', 'CITED2', 'CLSPN', 'GINS2', 'CCNE2', 'H4C3', 'H1-2', 'H1-4', 'PTTG1', 'C21orf58', 'CENPA', 'TUBA1B', 'STMN1', 'BIRC5')
# convert to ensembl IDs
ensembl_important = mapIds(org.Hs.eg.db, keys = important, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
#ensembl_features3 = na.omit(data.frame(ensembl_features3))
#ensembl_features3_plot = ensembl_features3$ensembl_features3

g1_important = c('HMGN2', 'SMS', 'HMGB1', 'IGFBP5', 'STMN1', 'DTYMK', 'CKS1B', 'CLU', 'NUCKS1', 'PRCP', 'MIAT', 'PTN', 'CHCHD2', 'GGCT', 'CCN2')
ensembl_g1 = mapIds(org.Hs.eg.db, keys = g1_important, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')

# Find marker genes
cluster_markers = read.csv(file.path(resdir2, 'analysis_output/', paste0(tag1,'_scTransform_Markers_together.csv')))
cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

# Subset to marker genes in ccAFv2 marker genes
top10[top10$X %in% mgenes,]$X

# Visualize UMAP
p1 = DimPlot(seurat1, reduction = "umap", label=F, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub2])
c1 = DimPlot(seurat1, reduction = "umap", label=F, pt.size = 1, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1])
lst = list(p1, c1)

seurat_subset = subset(seurat1, subset = ccAFv2 != 'Unknown')

#pdf(file.path(resdir2, paste0('W9-1_heatmap_SCT.pdf')), height = 8, width = 10)
pdf(file.path(resdir2, paste0('W9-1_heatmap_SCT_important_features.pdf')), height = 12, width = 10)
#DoHeatmap(subset(seurat_subset,downsample = 1000), group.by = 'ccAFv2', features=goi_lst_ensembl2, size = 3, group.colors  = ccAFv2_colors) + scale_y_discrete(labels = rev(names(goi_lst_ensembl2)))
#DoHeatmap(seurat_subset, group.by = 'ccAFv2', features=top10[top10$X %in% mgenes,]$X, size = 3, group.colors  = ccAFv2_colors) + scale_y_discrete(labels = rev(top10[top10$X %in% mgenes,]$gene))
DoHeatmap(seurat_subset, group.by = 'ccAFv2', features=ensembl_important, group.colors  = ccAFv2_colors) + scale_y_discrete(labels = rev(names(ensembl_important)))
dev.off()

#pdf(file.path(resdir2, paste0('W9-1_dotplot_SCT.pdf')), height = 15, width = 45)
pdf(file.path(resdir2, paste0('W9-1_dotplot_SCT_G1_important_features.pdf')), height = 15, width = 45)
DotPlot(seurat1, group.by = 'ccAFv2', features=ensembl_g1) + RotatedAxis()
dev.off()

#AverageExpression(seurat_subset, features = ensembl_features3_plot, group.by = 'ccAFv2')


# Read in other calls
tricycle_calls = read.csv('compare_classifiers/PCW9/PCW9_tricycle_calls.csv', row.names = 'X')
schwabe_calls = read.csv('compare_classifiers/PCW9/PCW9_ccschwabe_calls.csv', row.names = 'X')
recat_calls = read.csv('compare_classifiers/PCW9/PCW9_recat_calls.csv', row.names = 'X')
peco_calls = read.csv('compare_classifiers/PCW9/PCW9_peco_calls.csv', row.names = 'X')
ccaf_calls = read.csv('compare_classifiers/PCW9/PCW9_ccaf_calls.csv', row.names = 'X')
cyclone_calls = read.csv('compare_classifiers/PCW9/PCW9_cyclone_calls.csv', row.names = 'X')

seurat1$tricycle = tricycle_calls$x
seurat1$ccschwabe = schwabe_calls$x
seurat1$recat = recat_calls$x
seurat1$peco = peco_calls$x
seurat1$ccaf = ccaf_calls$ccAF
seurat1$cyclone = cyclone_calls$x

# ccAF
sub3 = ccAF_order %in% factor(seurat1$ccaf)
seurat1$ccaf <- factor(seurat1$ccaf, levels = ccAF_order[sub3])
# Tricycle
sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 = factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])
sub2 = ccSeurat_order %in% factor(seurat1$Phase)
seurat1$Phase <- factor(seurat1$Phase, levels = ccSeurat_order[sub2])
sub3 = ccAF_order %in% factor(seurat1$ccaf)
seurat1$ccaf = factor(seurat1$ccaf, levels = ccAF_order[sub3])
sub4 = tricycle_order %in% factor(seurat1$tricycle)
seurat1$tricycle <- factor(seurat1$tricycle, levels = tricycle_order[sub4])
# ccSchwabe
sub5 = ccSchwabe_order %in% factor(seurat1$ccschwabe)
seurat1$ccschwabe <- factor(seurat1$ccschwabe, levels = ccSchwabe_order[sub5])
# reCAT
sub6 = recat_order %in% factor(seurat1$recat)
seurat1$recat <- factor(seurat1$recat, levels = recat_order[sub6])
# cyclone
# same as Seurat
sub7 = ccSeurat_order %in% factor(seurat1$cyclone)
seurat1$cyclone <- factor(seurat1$cyclone, levels = ccSeurat_order[sub7])
# peco
# same as tricycle
sub8 = tricycle_order %in% factor(seurat1$peco)
seurat1$peco <- factor(seurat1$peco, levels = tricycle_order[sub8])


# Downstream analysis
seurat1 = RunPCA(seurat1, dims = 1:30, verbose=FALSE)
seurat1 = FindNeighbors(seurat1, dims = 1:30, verbose=FALSE)
seurat1 = FindClusters(seurat1, verbose=FALSE, resolution = 0.8)
seurat1 = RunUMAP(seurat1, dims=1:30, verbose=FALSE)

#--------------------------------
# Plot
#--------------------------------

# Visualize umap
d1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6,  pt.size = 1)
c1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccaf', cols = ccAF_colors[sub3])
c2 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccAFv2', cols = ccAFv2_colors[sub1])
p1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub2])
t1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'tricycle', cols = tricycle_colors[sub4])
s1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccschwabe', cols = ccSchwabe_colors[sub5])
r1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'recat', cols = recat_colors[sub6])
cy1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'cyclone', cols = ccSeurat_colors[sub3])
pe1 = DimPlot(seurat1, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'peco', cols = tricycle_colors[sub4])

# Plot
tag = 'PCW9'
pdf(file.path('compare_classifiers/PCW9', paste0(tag, '_classifiers_umap_visualization.pdf')), height = 10, width = 24)
#pdf(file.path(resdir2, paste0(tag, '_classifiers_umap_visualization.pdf')), height = 12, width = 24)
lst1 = list(c1, c2, p1, t1, r1, s1, pe1, cy1)
#lst1 = list(c2, c1, p1, t1, r1, s1, pe1, d1)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2, 3, 4), c(5,6,7, 8)), top = "")
dev.off()
