##########################################################
## ccAFv2:  Plotting cell cycle classifiers U5          ##
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

#--------------------------------
# Set up section / load packages
#--------------------------------
#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccafv2_extra

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)

#setwd("files/")

# Plotting order & colors
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccSeurat_order = c('G1', 'S', 'G2M')
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")
tricycle_order = c('G1/G0', 'S', 'G2/M', 'M')
tricycle_colors = c("G1/G0" = "#f37f73", "S" = "#8571b2", "G2/M" = "#3db270", "M" = "#6d90ca", "Unknown" = "#D3D3D3")
ccSchwabe_order = c('G1.S', 'S', 'G2', 'G2.M', 'M.G1', 'Unknown')
ccSchwabe_colors = c("G1.S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2.M" = "#3db270", "M.G1"= "#6d90ca", "Unknown" = "#D3D3D3")
recat_order = c('G1', 'G1S', 'S', 'G2', 'G2M', 'M')
recat_colors = c("G1" = "#f37f73", "G1S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2M" = "#3db270", "M"= "#6d90ca")

# Set directory structures
tag = 'U5'
resdir = 'data'
savedir = 'compare_classifiers'
data_dir = file.path(resdir, tag)
resdir2 = file.path(savedir, tag)

# Load in data
seurat2 = readRDS(file.path(data_dir, paste0(tag, '_normalized_gene_symbols.rds')))

# Load in all classifier calls
classifiers = c('ccafv2', 'seurat', 'tricycle', 'ccschwabe', 'recat', 'cyclone', 'peco')
for (class1 in classifiers){
  calls = read.csv(file.path(resdir2, paste0(tag, '_', class1, '_calls.csv')), row.names = 1)
  colnames(calls) = c('x')
  seurat2[[class1]] = calls$x
}

features <- c("CCND1")

r1 = RidgePlot(tmp2, features, cols = c("#1fb1a9"))
r2 = RidgePlot(tmp3, features, cols = c("#f37f73"))
r3 = RidgePlot(tmp4, features, cols = c("#1fb1a9"))
lst1 = list(r2, r1, r3)

pdf('ccAFv2_Late_G1_ccAF_G1_cells.pdf', width = 16, height = 4)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2, 3)), top = "")
DotPlot(tmp2, features = c('CCND1'))
DotPlot(tmp3, features = c('CCND1'))
DotPlot(tmp4, features = c('CCND1', 'CCNE2'))
dev.off()



tmp3 = subset(seurat2, subset = ccafv2 == 'G1')
tmp4 = subset(seurat2, subset = ccafv2 == 'Late G1')

#--------------------------------
# Organize states for each classifer
#--------------------------------
# ccAF & ccAFv2
sub1 = ccAFv2_order %in% factor(seurat2$ccAF)
seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub1])
sub2 = ccAFv2_order %in% factor(seurat2$ccafv2)
seurat2$ccafv2 <- factor(seurat2$ccafv2, levels = ccAFv2_order[sub2])
# Seurat
sub3 = ccSeurat_order %in% factor(seurat2$seurat)
seurat2$seurat <- factor(seurat2$seurat, levels = ccSeurat_order[sub3])
# Tricycle
sub4 = tricycle_order %in% factor(seurat2$tricycle)
seurat2$tricycle <- factor(seurat2$tricycle, levels = tricycle_order[sub4])
# ccSchwabe
sub5 = ccSchwabe_order %in% factor(seurat2$ccschwabe)
seurat2$ccschwabe <- factor(seurat2$ccschwabe, levels = ccSchwabe_order[sub5])
# reCAT
sub6 = recat_order %in% factor(seurat2$recat)
seurat2$recat <- factor(seurat2$recat, levels = recat_order[sub6])
# cyclone
# same as Seurat
sub7 = ccSeurat_order %in% factor(seurat2$cyclone)
seurat2$cyclone <- factor(seurat2$cyclone, levels = ccSeurat_order[sub7])
# peco
# same as tricycle
sub8 = tricycle_order %in% factor(seurat2$peco)
seurat2$peco <- factor(seurat2$peco, levels = tricycle_order[sub8])

#--------------------------------
# Plot
#--------------------------------
# Visualize umap
d1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6,  pt.size = 1)
c1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccAF', cols = ccAFv2_colors[sub1])
c2 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccafv2', cols = ccAFv2_colors[sub2])
p1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'seurat', cols = ccSeurat_colors[sub3]) + ggtitle('ccSeurat')
t1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'tricycle', cols = tricycle_colors[sub4])
s1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccschwabe', cols = ccSchwabe_colors[sub5])
r1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'recat', cols = recat_colors[sub6])
cy1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'cyclone', cols = ccSeurat_colors[sub3])
pe1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'peco', cols = tricycle_colors[sub4])

# Plot
pdf(file.path(resdir2, paste0(tag, '_classifiers_umap_visualization.pdf')), height = 12, width = 24)
lst1 = list(c1, c2, p1, t1, r1, s1, pe1, cy1)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2, 3, 4), c(5,6,7, 8)), top = "")
dev.off()
