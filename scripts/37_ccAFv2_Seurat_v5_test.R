##########################################################
## ccAFv2:  testing ccAFv2 with seurat v5               ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccafv2_seurat5

remotes::install_version("matrixStats", version="1.1.0")
library(matrixStats)
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

# Set working directory
setwd("files/")
tag1 = 'NSC'
resdir = 'testData/GSE155121'
savedir = file.path(resdir, paste0(tag1))


# Plotting order & colors
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

# Load NSC data
seurat1 <- LoadH5Seurat(file.path(resdir,'nsc.h5seurat'))
seurat1
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito = colSums(seurat1[mito.genes,])/colSums(seurat1@assays[['RNA']])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")
nCount = colSums(x = seurat1, slot = "counts")  # nCount_RNA
seurat1 <- AddMetaData(seurat1, nCount, col.name = "nCount_RNA")

# Subset to W8-1 week stage
week_stage = 'W8-1'
h1 = 0.001
h2 = 0.04
v1 = 4000
v2 = 23000

ws1 = week_stage
cat('\n',week_stage,'\n')
seurat_subset = subset(seurat1, subset = week_stage == ws1)
cat('\n', dim(seurat_subset),'\n')
# Quality control filtering
keep.detect = which(seurat_subset@meta.data$percent.mito < h2 & seurat_subset@meta.data$percent.mito > h1 & seurat_subset@meta.data$nCount_RNA < v2 & seurat_subset@meta.data$nCount_RNA > v1)
seurat_subset = subset(seurat_subset, cells=colnames(seurat_subset)[keep.detect])
cat('\n', dim(seurat_subset),'\n')
cat('Normalization \n')
# Change to Ensembl ID
seurat_subset@assays$RNA@data@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
seurat_subset@assays$RNA@counts@Dimnames[1][[1]] = seurat_subset@assays$RNA@meta.features$gene_ids
seurat_subset = SCTransform(seurat_subset, return.only.var.genes = FALSE, vst.flavor = 'v1') # SCTransform v4
seurat_subset = PredictCellCycle(seurat_subset, do_sctransform=FALSE)
seurat_subset$ccAFv2 = as.character(seurat_subset$ccAFv2)
seurat_subset$ccAFv2_v5 = seurat_subset$ccAFv2


# Read in ccAFv2_R and ccAFv2_py predictions on same data
preds = read.csv('R_vs_python_GSE155121_compare.csv', row.names = 'X')
df = data.frame(preds$ccAFv2, preds$ccAFv2_py, seurat_subset$ccAFv2_v5)
colnames(df) = c('ccAFv2_R', 'ccAFv2_py', 'ccAFv2_seurat_v5')

# Get ccAFv2 frequencies for each version
tmp = data.frame(table(df$ccAFv2_R))
tmp2 = data.frame(table(df$ccAFv2_py))
tmp3 = data.frame(table(df$ccAFv2_seurat_v5))
tmp$version = 'ccAFv2_R'
tmp2$version = 'ccAFv2_py'
tmp3$version = 'ccAFv2_seurat_v5'

# Bind all together and plot
df2 = bind_rows(tmp, tmp2, tmp3)
df2$Var1 = factor(df2$Var, levels = ccAFv2_order)
pdf(file.path('ccAFv2_R_python_v5.pdf'))
ggplot(df2, aes(fill=Var1, y=Freq, x=version)) + scale_fill_manual(values=ccAFv2_colors) +
    geom_bar(position="fill", stat="identity")
dev.off()
