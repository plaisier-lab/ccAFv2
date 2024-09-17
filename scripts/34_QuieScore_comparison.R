##########################################################
## ccAFv2:  QuieScore comparison                        ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccafv2_seurat4

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
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(pheatmap)

library(devtools)
install_github("https://github.com/dkornai/QuieScore")
library(QuieScore)

setwd("files/")

# Plotting order & colors
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

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

# Prepare for QuieScore
mat.expr = data.frame(seurat2[['RNA']]@data)
processedData <- processInput(mat.expr, cancer_type = "LGG", gene_naming = "name", log_transformed=FALSE)
#processedData <- processInput(mat.expr, cancer_type = "GBM", gene_naming = "name", log_transformed=FALSE)

G0scores <- QuiescenceScore(processedData)
G0scores$G0 <- ifelse(G0scores$q_score_raw>3, 'Neural G0', 'NA')
seurat2$G0 = G0scores$G0

# Order ccAFv2 calls
sub1 = ccAFv2_order %in% factor(seurat2$ccAF)
seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub1])
sub2 = ccAFv2_order %in% factor(seurat2$ccafv2)
seurat2$ccafv2 <- factor(seurat2$ccafv2, levels = ccAFv2_order[sub2])

# Functions
plotHyper = function(data, tag, save_dir = file.path(resdir2)){
  datas = data
  savedir = save_dir
  compare_col = group.by
  df = table(datas$ccafv2, datas$G0)
  colnames(df) = c('not G0', 'G0')
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
  #pm[pm>20] = 20
  pdf(file.path(savedir, paste0(tag, '_ccAFv2_quiescore_heatmap_091124.pdf')), height = 8, width = 8)
  print(pheatmap(pm, cluster_cols = F, cluster_rows = F, colorRampPalette(c("white", "red"))(100), display_numbers = round(-log10(as.matrix(enrich)),2)))
  dev.off()
}

# Plot
plotHyper(seurat2, 'U5')
