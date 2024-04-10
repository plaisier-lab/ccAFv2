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
#sessionInfo()
library(ssgsea.GBM.classification)
#library(verification)
#library(MCMCpack)

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

# Set working directory
setwd("files/")
tag = 'U5'
resdir2 = file.path('U5/analysis_output')
resdir4 = file.path('U5/analysis_output/seurat_objects/final')


#--------------------------------
# Load in data
#--------------------------------
# Load in data
seurat2 = readRDS(file.path(resdir4, paste0(tag, "_normalized_gene_symbols.rds")))

# Read in 600, 200 v2 model ccAFv2 calls and add as metadata
ccAFv2_calls = read.csv('U5_ccAFv2_calls_010424.csv', row.names = 1)
seurat2$ccAFv2 = ccAFv2_calls$x

# Read in Phase calls and add as metadata
ccseurat <- read_csv('U5_ccseurat_calls_101023.csv')
seurat2$Phase = ccseurat$x

# Read in tricycle calls and add as metadata
tricycle_calls = read.csv('U5_tricycle_calls_020224.csv', row.names = 1)
seurat2$tricycle = tricycle_calls$x

# Read in ccSchwabe calls and add as metadata
schwabe <- read.csv('U5_schwabe_calls_101023.csv', row.names = 'X')
seurat2$ccSchwabe = schwabe$x

# Read in peco calls and add as metadata
peco_calls = read.csv('U5_peco_calls_020224.csv', row.names = 1)
seurat2$peco = peco_calls$x

# Read in recat calls and add as metadata
recat <- read_csv('U5_recat_calls_101023.csv')
seurat2$recat = recat$recat_calls

# Read in cyclone calls and add as metadata
cyclone <- read_csv('U5_cyclone_calls_101023.csv')
seurat2$cyclone = cyclone$results.phases


#--------------------------------
# Organize states for each classifer
#--------------------------------
# ccAF & ccAFv2
sub1 = ccAFv2_order %in% factor(seurat2$ccAF)
seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub1])
sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
# Seurat
sub3 = ccSeurat_order %in% factor(seurat2$Phase)
seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub3])
# Tricycle
sub4 = tricycle_order %in% factor(seurat2$tricycle)
seurat2$tricycle <- factor(seurat2$tricycle, levels = tricycle_order[sub4])
# ccSchwabe
sub5 = ccSchwabe_order %in% factor(seurat2$ccSchwabe)
seurat2$ccSchwabe <- factor(seurat2$ccSchwabe, levels = ccSchwabe_order[sub5])
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
c2 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccAFv2', cols = ccAFv2_colors[sub2])
p1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub3]) + ggtitle('ccSeurat')
t1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'tricycle', cols = tricycle_colors[sub4])
s1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'ccSchwabe', cols = ccSchwabe_colors[sub5])
r1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'recat', cols = recat_colors[sub6])
cy1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'cyclone', cols = ccSeurat_colors[sub3])
pe1 = DimPlot(seurat2, reduction = "umap", label=F, label.size = 6, pt.size = 1, group.by = 'peco', cols = tricycle_colors[sub4])

# Plot
pdf(file.path(resdir2, paste0(tag, "_visualization_test_020224.pdf")), height = 12, width = 24)
lst1 = list(c1, c2, p1, t1, r1, s1, pe1, cy1)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2, 3, 4), c(5,6,7, 8)), top = "")
dev.off()
