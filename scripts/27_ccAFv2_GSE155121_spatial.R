##########################################################
## ccAFv2:  testing ccAFv2 on ST                        ##
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


# docker run -it -v '/home/soconnor/old_home/ccNN/testData:/files' cplaisier/ccafv2_extra

library(reticulate)
use_python('/usr/bin/python3')

#devtools::install_github('satijalab/seurat-data')
devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ccAFv2)

# Plotting order & colors
ccAF_colors = c("G1/other" = "#9aca3c", "Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
ccAF_order = c("G1/other", 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccSeurat_order = c('G1', 'S', 'G2M')

setwd('/files/GSE155121_spatial')

samples = c("GSM6736780_Spatial_10x_PCW4_20220122_slice1.RData", "GSM6736781_Spatial_10x_PCW4_20220122_slice2.RData", "GSM6736782_Spatial_10x_PCW4_20220122_slice3.RData", "GSM6736783_Spatial_10x_PCW4_20220122_slice4.RData", "GSM6736784_Spatial_10x_PCW4_20220122_slice5.RData", "GSM6736785_Spatial_10x_PCW4_20220122_slice6.RData", "GSM6736786_Spatial_10x_PCW4_20220122_slice7.RData", "GSM6736787_Spatial_10x_PCW4_20220122_slice8.RData")

for(sample in samples) {
    load(sample)
    name1 = gsub('.RData', '', sample)
    GW7 = PredictCellCycle(GW6, species='human', gene_id='symbol', spatial=TRUE)
    saveRDS(GW7, file=paste0('results/',gsub('RData','rds', sample)))
    pdf(paste0('results/',name1,'_wCellCycle.pdf'))
    s1 = SpatialFeaturePlot(GW7, features=NULL, alpha=0)
    s2 = SpatialFeaturePlot(GW7, features = "nCount_Spatial") + theme(legend.position = "right")
    s3 = SpatialFeaturePlot(GW7, features=c('MKI67','MCM2'))
    s4 = SpatialDimPlot.ccAFv2(GW7) + theme(legend.position = "right")
    print(s1)
    print(s2)
    print(s3)
    print(s4)
    dev.off()
}

samples = c("GSM6736780_Spatial_10x_PCW4_20220122_slice1", "GSM6736781_Spatial_10x_PCW4_20220122_slice2", "GSM6736782_Spatial_10x_PCW4_20220122_slice3", "GSM6736783_Spatial_10x_PCW4_20220122_slice4", "GSM6736784_Spatial_10x_PCW4_20220122_slice5", "GSM6736785_Spatial_10x_PCW4_20220122_slice6", "GSM6736786_Spatial_10x_PCW4_20220122_slice7", "GSM6736787_Spatial_10x_PCW4_20220122_slice8")
datas = list()
for(sample in samples){
  datas[[sample]] = readRDS(file.path(paste0('results/',sample, '.rds')))
  slice_tag = unlist(strsplit(sample, '_'))[6]
  datas[[sample]][['slice']] = slice_tag
  print(length(unique(intersect(datas[[sample]]@assays$SCT@data@Dimnames[[1]], mgenes))))
}

# merge
merged = merge(x = datas[['GSM6736780_Spatial_10x_PCW4_20220122_slice1']], y = c(datas[['GSM6736781_Spatial_10x_PCW4_20220122_slice2']], datas[['GSM6736782_Spatial_10x_PCW4_20220122_slice3']], datas[['GSM6736783_Spatial_10x_PCW4_20220122_slice4']], datas[['GSM6736784_Spatial_10x_PCW4_20220122_slice5']], datas[['GSM6736785_Spatial_10x_PCW4_20220122_slice6']], datas[['GSM6736786_Spatial_10x_PCW4_20220122_slice7']], datas[['GSM6736787_Spatial_10x_PCW4_20220122_slice8']]))
dim(merged) #[1] 19499 14284

seurat1 <- merged
subset1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 <- factor(seurat1$ccAFv2, levels = ccAFv2_order[subset1])

#--- ccAFv2 vs. week stage stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$slice)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cf_1 = matrix(ncol=length(unique(seurat1$slice)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$slice)))){
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cf)
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf(paste0('spatial_ccAFv2_percentages_020624.pdf'))
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 15, y=max(colSums(cf_1)), bty = "n"))
dev.off()
