##########################################################
## ccAFv2:  testing ccAFv2                              ##
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

### Docker command to start-up Seurat capable analysis container
# docker run -it -v '/home/soconnor/old_home:/files' cplaisier/ccafv2_extra

### Packages required to run analyses
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

setwd('/files')

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)

# Plotting order & colors
ccAF_colors = c("G1/other" = "#9aca3c", "Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
ccAF_order = c("G1/other", 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccSeurat_colors = c("G1"="#f37f73", "S"="#8571b2", "G2M"="#3db270")
ccSeurat_order = c('G1', 'S', 'G2M')

###########################
### Nowakowski_2017 GBM ###
###########################

# Load up count data and start Seurat
setwd('/files/scGlioma/Nowakowski')
d1 = read.table(gzfile('exprMatrix.tsv.gz'),sep='\t',row.names=1,header=T)
dim(d1) #[1] 56864  4261
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1) #[1] 49870  4261
d1 = na.omit(d1)

# Make Seurat 3 object
Nowakowski_2017 = CreateSeuratObject(counts=d1, min.cells = 3, min.features = 200)
dim(Nowakowski_2017) #[1] 41222  4261

# Load meta-data
meta = read.table('meta.tsv', header=T, row.names=1, sep='\t', as.is=T)
meta[is.na(meta)] = 'none'

# Make a cell type conversion list
conversion = list()
conversion[['Astrocyte']] = 'Astrocyte'
conversion[['Choroid']] = 'Choroid'
conversion[['Endothelial']] = 'Endothelial'
conversion[['EN-PFC1']] = 'EN-PFC'
conversion[['EN-PFC2']] = 'EN-PFC'
conversion[['EN-PFC3']] = 'EN-PFC'
conversion[['EN-V1-1']] = 'EN-V1'
conversion[['EN-V1-2']] = 'EN-V1'
conversion[['EN-V1-3']] = 'EN-V1'
conversion[['Glyc']] = 'Glyc'
conversion[['IN-CTX-CGE1']] = 'IN-CTX-CGE'
conversion[['IN-CTX-CGE2']] = 'IN-CTX-CGE'
conversion[['IN-CTX-MGE1']] = 'IN-CTX-MGE'
conversion[['IN-CTX-MGE2']] = 'IN-CTX-MGE'
conversion[['IN-STR']] = 'IN-STR'
conversion[['IPC-div1']] = 'IPC-div'
conversion[['IPC-div2']] = 'IPC-div'
conversion[['IPC-nEN1']] = 'IPC-nEN'
conversion[['IPC-nEN2']] = 'IPC-nEN'
conversion[['IPC-nEN3']] = 'IPC-nEN'
conversion[['MGE-div']] = 'MGE-div'
conversion[['MGE-IPC1']] = 'MGE-IPC'
conversion[['MGE-IPC2']] = 'MGE-IPC'
conversion[['MGE-IPC3']] = 'MGE-IPC'
conversion[['MGE-RG1']] = 'MGE-RG'
conversion[['MGE-RG2']] = 'MGE-RG'
conversion[['Microglia']] = 'Microglia'
conversion[['Mural']] = 'Mural'
conversion[['none']] = 'none'
conversion[['nEN-early1']] = 'nEN-early'
conversion[['nEN-early2']] = 'nEN-early'
conversion[['nEN-late']] = 'nEN-late'
conversion[['nIN1']] = 'nIN'
conversion[['nIN2']] = 'nIN'
conversion[['nIN3']] = 'nIN'
conversion[['nIN4']] = 'nIN'
conversion[['nIN5']] = 'nIN'
conversion[['OPC']] = 'OPC'
conversion[['oRG']] = 'oRG'
conversion[['RG-div1']] = 'RG-div'
conversion[['RG-div2']] = 'RG-div'
conversion[['RG-early']] = 'RG-early'
conversion[['tRG']] = 'tRG'
conversion[['U1']] = 'U'
conversion[['U2']] = 'U'
conversion[['U3']] = 'U'
conversion[['U4']] = 'U'
conversion[['vRG']] = 'vRG'
meta[,'WGCNAcluster_restricted'] = unlist(sapply(meta[,'WGCNAcluster'], function(x) { ifelse(is.na(x), 'NA', conversion[[x]]) }))

# Add meta data variables of interest
for(pheno1 in colnames(meta)) {
    Nowakowski_2017[[pheno1]] = meta[,pheno1]
}

# Classify with ccAFv2
seurat1 = Nowakowski_2017
seurat1 = PredictCellCycle(seurat1, assay = 'RNA', gene_id = 'symbol')
table(seurat1$ccAFv2)

write.csv(seurat1$ccAFv2, file.path('Nowakowski_ccAFv2_calls.csv'))
write.csv(table(seurat1$ccAFv2, seurat1$WGCNAcluster_restricted), file.path('Nowakowski_ccAFv2_calls_by_cell_type.csv'))
saveRDS(seurat1, file.path('Nowakowski_2017_with_ccAFv2_calls.rds'))

# Prepare for plotting
df1 = data.frame(table(seurat1$ccAFv2))
rownames(df1) = df1$Var1
df1$perc = df1$Freq/dim(seurat1)[2]*100
sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
seurat1$ccAFv2 <- factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])
# organize x-axis order
seurat1$WGCNAcluster_restricted <- factor(seurat1$WGCNAcluster_restricted, levels = c('Astrocyte', 'OPC', 'oRG', 'tRG', 'vRG', 'MGE-RG', 'RG-early', 'RG-div', 'Choroid', 'Microglia', 'IN-STR', 'IN-CTX-MGE', 'IN-CTX-CGE', 'nIN', 'EN-PFC', 'EN-V1', 'nEN-late', 'nEN-early', 'IPC-nEN', 'IPC-div', 'MGE-IPC', 'MGE-div', 'Endothelial', 'Mural', 'Glyc','U', 'none'))

#--- ccAFv2 vs. cluster ids stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$WGCNAcluster_restricted)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cf_1 = matrix(ncol=length(unique(seurat1$WGCNAcluster_restricted)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$WGCNAcluster_restricted)))){
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cf)
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf('ccAFv2_percentages_010824.pdf', width = 10, height = 8)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 15, y=max(colSums(cf_1)), bty = "n"))
dev.off()


#################################
## GSE67833: Llorens-Bobadilla ##
#################################

setwd('/files/scGlioma/GSE67833')
exp_mat = read.csv('GSE67833_Gene_expression_matrix.csv',header=T,row.names=1)

GSE67833 = CreateSeuratObject(counts = exp_mat, min.cells = 3, min.features = 200)
dim(GSE67833)
GSE67833 = NormalizeData(object = GSE67833)
GSE67833 = FindVariableFeatures(object = GSE67833, selection.method='vst', nfeatures=3000)
GSE67833 = ScaleData(GSE67833)

meta_data = read.csv('meta_data.csv', header=T, row.names=2)
GSE67833[['cell.type']] = meta_data[colnames(exp_mat), 2]

# For classification
GSE67833.NI = subset(GSE67833, cells = colnames(GSE67833)[-grep('Injured',GSE67833[['cell.type']][,1])])
GSE67833.NI = subset(GSE67833, cells = colnames(GSE67833.NI)[!is.na(GSE67833.NI[['cell.type']])])

cutoff1 = 0.5
# Classify with ccAFv2
GSE67833.NI = PredictCellCycle(GSE67833.NI, assay = 'RNA', cutoff = cutoff1, species = 'mouse', gene_id='ensembl')
saveRDS(GSE67833.NI, file.path('Llorens_Bobadilla_with_ccAFv2_calls.rds'))

# Prepare for plotting
seurat1 = GSE67833.NI
df1 = data.frame(table(seurat1$ccAFv2))
rownames(df1) = df1$Var1
df1$perc = df1$Freq/dim(seurat1)[2]*100

sub1 = ccAFv2_order %in% factor(seurat1$ccAFv2)
# organize ccAFv2
seurat1$ccAFv2 <- factor(seurat1$ccAFv2, levels = ccAFv2_order[sub1])
# organize x-axis order
seurat1$cell.type <- factor(seurat1$cell.type, levels = c('qNSC1', 'qNSC2', 'aNSC1', 'aNSC2', 'Olig', 'NB'))

#--- ccAFv2 vs. cluster ids stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$cell.type)
cf <- cf[rowSums(cf[])>0,]
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cf_1 = matrix(ncol=length(unique(seurat1$cell.type)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$cell.type)))){
  print(i)
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cf)
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf(paste0('GSE67833_ccAFv2_percentages_011224_', cutoff1,'.pdf'), width = 10, height = 8)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 15, y=max(colSums(cf_1)), bty = "n"))
dev.off()


############################
### Dulken: PRJNA324289 ####
############################
setwd('/files/scGlioma/PRJNA324289')

# Load in data
d1 = read.csv('Counts_AllLiveSingleCells_IN_VIVO_ONLY.csv',header=T,row.names=1)
ct1 = sapply(colnames(d1), function(x) { strsplit(x,'_')[[1]][1] })
PRJNA324289 = CreateSeuratObject(counts = d1, min.cells = 3, min.features = 200)
PRJNA324289[['cell_type']] = ct1

# Classify with ccAFv2
cutoff1 = 0.5
PRJNA324289 = PredictCellCycle(PRJNA324289, species='mouse', gene_id='symbol', cutoff = cutoff1)
saveRDS(PRJNA324289, file.path('Dulken_with_ccAFv2_calls.rds'))

# Prepare for plotting
#df1 = data.frame(table(PRJNA324289$ccAFv2))
#rownames(df1) = df1$Var1
#df1$Freq/dim(PRJNA324289)[2]*1000
#df1$perc = df1$Freq/dim(PRJNA324289)[2]*100

seurat1 = PRJNA324289
# organize x-axis order
seurat1$cell_type <- factor(seurat1$cell_type, levels = c('qNSC', 'aNSC', 'Ast', 'NPC'))

#--- ccAFv2 vs. cluster ids stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$cell_type)
cf <- cf[rowSums(cf[])>0,]
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cf_1 = matrix(ncol=length(unique(seurat1$cell_type)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$cell_type)))){
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cf)
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf(paste0('Dulken_ccAFv2_percentages_011024_', cutoff1,'.pdf'), width = 10, height = 8)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 15, y=max(colSums(cf_1)), bty = "n"))
dev.off()


#################################
##  GSE165555 ##
#################################

setwd('/files/ccNN/testData/GSE165555')
# Load scRNA-seq data
data <- readRDS('GSM5039270_scSeq.rds')
dim(data) #25803 24261
# Classify with ccAFv2
data = PredictCellCycle(data, species='mouse', gene_id='symbol')
saveRDS(data, file.path('Cebrian_Silla_with_ccAFv2_calls.rds'))
seurat1 = data
df1 = data.frame(table(seurat1$ccAFv2))
rownames(df1) = df1$Var1
df1$perc = df1$Freq/dim(seurat1)[2]*100


# Prepare for plotting
# organize x-axis order
seurat1$Cell_Type <- factor(seurat1$Cell_Type, levels = c('Ependymal cells', 'Astrocytes', 'B cells', 'C cells', 'A cells', 'GABAergic neurons', 'Neuron', 'OPC/Oligo', 'Microglia', 'Endothelial cells', 'Pericytes/VSMC', 'VLMC1', 'Mitosis'))
#--- ccAFv2 vs. cluster ids stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$Cell_Type)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cnewdf = cnewdf[apply(cnewdf[,-1], 1, function(x) !all(x==0)),]
cf_1 = matrix(ncol=length(unique(seurat1$Cell_Type)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$Cell_Type)))){
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cnewdf)[-length(rownames(cnewdf))]
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf('GSE165555_ccAFv2_percentages_010824.pdf', width = 10, height = 6)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 8, y=max(colSums(cf_1)), bty = "n"))
dev.off()

#--- ccSeurat vs. cluster ids stacked barplot ---#
cf <- table(seurat1$Phase, seurat1$Cell_Type)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cnewdf = cnewdf[apply(cnewdf[,-1], 1, function(x) !all(x==0)),]
cf_1 = matrix(ncol=length(unique(seurat1$Cell_Type)), nrow=length(unique(seurat1$Phase)))
for(i in c(1:length(unique(seurat1$Cell_Type)))){
  for(n in c(1:length(unique(seurat1$Phase)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$Phase))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cnewdf)[-length(rownames(cnewdf))]
sub1 = rownames(data.frame(ccSeurat_colors)) %in% rownames(cf_1)
pdf('GSE165555_ccSeurat_percentages_010824.pdf', width = 10, height = 6)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccSeurat_colors[sub1], args.legend=list(x=ncol(cf_1) + 12, y=max(colSums(cf_1)), bty = "n"))
dev.off()

# Subset to just B cells & split by Predicted_Region_Bcells
seurat2 = subset(seurat1, subset = Cell_Type == 'B cells')
dim(seurat2) #[1] 25803  3369
seurat3 = subset(seurat2, subset = Predicted_Region_Bcells != 'NA')
dim(seurat3) #[1] 25803  2413
df = data.frame(seurat3$ccAFv2, seurat3$Predicted_Region_Bcells)
colnames(df) = c('ccAFv2', 'Predicted_Region')

pdf('GSE165555_B_cells_ccAFv2_by_region_010824.pdf', width = 10, height = 6)
ggplot(df, aes(ccAFv2)) + geom_bar(aes(fill = Predicted_Region), position = "dodge")
dev.off()


# Load snRNA-seq data
data1 <- readRDS('GSE165551_sNucSeq.rds')
dim(data1)# 23079 45820
# Classify with ccAFv2
data1 = PredictCellCycle(data1, species='mouse', gene_id='symbol')
seurat2 = data1
df1 = data.frame(table(seurat2$ccAFv2))
rownames(df1) = df1$Var1
df1$perc = df1$Freq/dim(seurat2)[2]*100

# organize x-axis order
seurat2$CellType_Detailed <- factor(seurat2$CellType_Detailed, levels = c('Ependymal', 'Astrocytes', 'B cells', 'Dividing/C cells', 'A cells', 'OPRM1 Neuron', 'VGLUT1 Neuron', 'CHAT Neuron', 'VGLUT2 Neuron', 'PTHLH Interneuron', 'PDYN Interneuron', 'Interneuron1', 'Interneuron2', 'SST Interneuron', 'TH Interneuron',  'Oligodendrocyte', 'OPC/Oligodendrocyte', 'Microglia', 'Endothelial', 'Pericyte/VSMC', 'SN-D1', 'SN-D2', 'Tac1 SN-D1', 'Vasc Leptomeningeal'))

#--- ccAFv2 vs. cluster ids stacked barplot ---#
cf <- table(seurat2$ccAFv2, seurat2$CellType_Detailed)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)
cnewdf = cnewdf[apply(cnewdf[,-1], 1, function(x) !all(x==0)),]
cf_1 = matrix(ncol=length(unique(seurat2$CellType_Detailed)), nrow=length(unique(seurat2$ccAFv2)))
for(i in c(1:length(unique(seurat2$CellType_Detailed)))){
  for(n in c(1:length(unique(seurat2$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cnewdf)[-length(rownames(cnewdf))]
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf('GSE165555_snrnaseq_ccAFv2_percentages_011824.pdf', width = 10, height = 6)
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 12, y=max(colSums(cf_1)), bty = "n"))
dev.off()


#################################
##  GSE136719.: Zhang 2021 ##
#################################
# Single‐cell analysis reveals dynamic changes of neural cells in developing human spinal cord

setwd('/files/ccNN/testData/GSE136719')

scProcessData = function(week1, location1, type1 = 'Cells', v1 = 3000, v2 = 50000, h1 = 0.001, h2 = 0.08, assay1 = 'SCT', cutoff1 = 0.5, resolution = 0.8, save_dir = 'analysis_output', obj_dir = 'seurat_objects', symbol = F){
  cat('\n Loading', week1, 'data\n')
  cat('\n Location:',location1,'\n')
  # Set directory to pull data from
  resdir1 = file.path(paste0(week1, '/', location1,'/',type1, '/outs/'))
  resdir2 = file.path(paste0(week1, '/', location1,'/', type1, '/', save_dir))
  resdir3 = file.path(paste0(week1, '/', location1, '/', type1,'/', obj_dir))
  # Make directories to save out analyses
  dir.create(resdir2, showWarnings = FALSE)
  dir.create(resdir3, showWarnings = FALSE)
  # Load data
  data = Read10X(file.path(resdir1), gene.column=2)
  rownames(data) = gsub("_", "-", rownames(data))
  cat('\n  Raw cells:', dim(data)[2], 'cells', dim(data)[1], 'genes \n')
  cat('\n Performing quality control\n')
  # Create seurat object
  seurat2 = CreateSeuratObject(counts = data, min.cells = 10, min.features = 200)
  seurat2$gw = week1
  cat('\n  Basic filter', dim(seurat2)[2], 'cells', dim(seurat2)[1], 'genes \n')
  mito_genes = grep('MT-', rownames(seurat2))
  seurat2[['percent.mito']] = PercentageFeatureSet(seurat2, features = mito_genes)/100
  # Plot QC
  pdf(file.path(resdir2, 'QC_plot_to_choose_cutoffs.pdf'))
  plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  # Filter based off QC plot
  keep.detect = which(seurat2@meta.data$percent.mito < h2 & seurat2@meta.data$percent.mito > h1 & seurat2@meta.data$nCount_RNA < v2 & seurat2@meta.data$nCount_RNA > v1)
  seurat2 = subset(seurat2, cells=colnames(seurat2)[keep.detect])
  cat('\n  Filtered to', dim(seurat2)[2], 'cells', dim(seurat2)[1], 'genes \n')
  # Check QC
  pdf(file.path(resdir2, 'QC_plot_to_choose_cutoffs_post.pdf'))
  plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  cat('\n Normalization\n')
  #seurat2 = SCTransform(seurat2, verbose = FALSE, return.only.var.genes = FALSE) # scale all genes
  #seurat2 = PredictCellCycle(seurat2, do_sctransform=FALSE, species='human', gene_id='symbol') # normalized by sct before so don't need to do it in PredictCellCycle function
  seurat2 = PredictCellCycle(seurat2, species='human', gene_id='symbol') #do normalization in function
  saveRDS(seurat2, file.path(paste0(resdir2, '/', week1, '_filtered_with_ccAFv2.rds')))
  seurat2 = SCTransform(seurat2, verbose = FALSE) # run regular sctransform after
  df1 = data.frame(table(seurat2$ccAFv2))
  df1 = df1[df1$Freq != 0,]
  rownames(df1) = df1$Var1
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% df1$Var1
  write.csv(data.frame(((df1['Freq']/dim(seurat2)[2])*100)), file.path(paste0(resdir2, '/', week1, '_ccAFv2_call_frequency.csv')))
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat2 = CellCycleScoring(seurat2, s.genes, g2m.genes, set.ident=FALSE)
  # Downstream analysis
  seurat2 <- RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 <- RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('\n Plotting UMAPs and ccAFv2 barplots \n')
  d1 = DimPlot(seurat2, reduction = "umap", label=F, group.by="seurat_clusters") + ggtitle("seurat_clusters")
  d2 = DimPlot(seurat2, reduction = "umap", label=F, group.by="Phase", cols = ccSeurat_colors) + ggtitle("Phase")
  d3 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAFv2", cols = ccAFv2_colors[sub1]) + ggtitle("ccAFv2")
  pdf(file.path(paste0(resdir2, '/', week1, '.pdf')), height = 8, width = 10)
  lst = list(d1, d2, d3)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = "")
  #--- ccAFv2 stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$gw)
  cf <- data.frame(cf[apply(cf!=0, 1, all),])
  colnames(cf) = c(week1)
  totals <- colSums(cf)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$gw)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$gw)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cnewdf)[-length(rownames(cnewdf))]
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(6, 9, 9, 6) + 2.0)
  print(barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1), col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 0.45, y=max(colSums(cf_1)), bty = "n")))
  dev.off()
  return(seurat2)
}

snProcessData = function(week1, location1, type1 = 'Nuc', v1 = 200, v2 = 15000, h1 = 200, h2 = 6000, assay1 = 'SCT', cutoff1 = 0.5, resolution = 0.8, save_dir = 'analysis_output', obj_dir = 'seurat_objects', symbol = F){
  cat('\n Loading', week1, 'data\n')
  cat('\n Location:',location1,'\n')
  # Set directory to pull data from
  resdir1 = file.path(paste0(week1, '/', location1,'/',type1, '/outs/'))
  resdir2 = file.path(paste0(week1, '/', location1,'/', type1, '/', save_dir))
  resdir3 = file.path(paste0(week1, '/', location1, '/', type1,'/', obj_dir))
  # Make directories to save out analyses
  dir.create(resdir2, showWarnings = FALSE)
  dir.create(resdir3, showWarnings = FALSE)
  # Load data
  data = Read10X(file.path(resdir1), gene.column=2)
  cat('\n  Raw cells:', dim(data)[2], 'cells', dim(data)[1], 'genes \n')
  rownames(data) = gsub("_", "-", rownames(data))
  cat('\n Performing quality control\n')
  # Create seurat object
  seurat2 = CreateSeuratObject(counts = data, min.cells = 10, min.features = 200)
  seurat2$gw = week1
  cat('\n  Basic filter', dim(seurat2)[2], 'cells', dim(seurat2)[1], 'genes \n')
  pdf(file.path(resdir2, 'QC_plot_to_choose_cutoffs.pdf'))
  plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$nFeature_RNA,
       xlab = 'nCount_RNA', ylab = 'nFeature_RNA', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
  dev.off()
  # Filter based off QC plot
  keep.detect = which(seurat2@meta.data$nFeature_RNA < h2 & seurat2@meta.data$nFeature_RNA > h1 & seurat2@meta.data$nCount_RNA < v2 & seurat2@meta.data$nCount_RNA > v1)
  seurat2 = subset(seurat2, cells=colnames(seurat2)[keep.detect])
  cat('\n  Filtered to', dim(seurat2)[2], 'cells', dim(seurat2)[1], 'genes \n')
  # Check QC
  pdf(file.path(resdir2, 'QC_plot_to_choose_cutoffs_post.pdf'))
  plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$nFeature_RNA,
       xlab = 'nCount_RNA', ylab = 'nFeature_RNA', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
  dev.off()
  cat('\n Normalization\n')
  #seurat2 = SCTransform(seurat2, verbose = FALSE, return.only.var.genes = FALSE) # scale all genes
  #seurat2 = PredictCellCycle(seurat2, do_sctransform=FALSE, species='human', gene_id='symbol') # normalized by sct before so don't need to do it in PredictCellCycle function
  seurat2 = PredictCellCycle(seurat2, species='human', gene_id='symbol') #do normalization in function
  saveRDS(seurat2, file.path(paste0(resdir2, '/', week1, '_filtered_with_ccAFv2.rds')))
  seurat2 = SCTransform(seurat2, verbose = FALSE) # run regular sctransform after
  df1 = data.frame(table(seurat2$ccAFv2))
  df1 = df1[df1$Freq != 0,]
  rownames(df1) = df1$Var1
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% df1$Var1
  write.csv(data.frame(((df1['Freq']/dim(seurat2)[2])*100)), file.path(paste0(resdir2, '/', week1, '_ccAFv2_call_frequency.csv')))
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat2 = CellCycleScoring(seurat2, s.genes, g2m.genes, set.ident=FALSE)
  # Downstream analysis
  seurat2 <- RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 <- FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 <- RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  #------------------------------------------------------
  # Plotting
  #---------------------------------------------------
  cat('\n Plotting UMAPs and ccAFv2 barplots \n')
  d1 = DimPlot(seurat2, reduction = "umap", label=F, group.by="seurat_clusters") + ggtitle("seurat_clusters")
  d2 = DimPlot(seurat2, reduction = "umap", label=F, group.by="Phase", cols = ccSeurat_colors) + ggtitle("Phase")
  d3 = DimPlot(seurat2, reduction = "umap", label=F, group.by="ccAFv2", cols = ccAFv2_colors[sub1]) + ggtitle("ccAFv2")
  pdf(file.path(paste0(resdir2, '/', week1, '.pdf')), height = 8, width = 10)
  lst = list(d1, d2, d3)
  grid.arrange(grobs = lst, layout_matrix = rbind(c(1, 2), c(3, NA)), top = "")
  #--- ccAFv2 stacked barplot ---#
  cf <- table(seurat2$ccAFv2, seurat2$gw)
  cf <- data.frame(cf[apply(cf!=0, 1, all),])
  colnames(cf) = c(week1)
  totals <- colSums(cf)
  cnewdf <- rbind(cf, totals)
  cf_1 = matrix(ncol=length(unique(seurat2$gw)), nrow=length(unique(seurat2$ccAFv2)))
  for(i in c(1:length(unique(seurat2$gw)))){
    for(n in c(1:length(unique(seurat2$ccAFv2)))) {
      cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat2$ccAFv2))+1, i]
    }
  }
  colnames(cf_1) = colnames(cf)
  rownames(cf_1) = rownames(cnewdf)[-length(rownames(cnewdf))]
  sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
  par(mar = c(6, 9, 9, 6) + 2.0)
  print(barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1), col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 0.45, y=max(colSums(cf_1)), bty = "n")))
  dev.off()
  return(seurat2)
}


datas_cells = list()
datas_cells[['GW8_S']] = scProcessData(week1 = 'GW8', location1 = 'Spinal_Whole')
datas_cells[['GW10_C']] = scProcessData(week1 = 'GW10', location1 = 'Cervical')
datas_cells[['GW10_L']] = scProcessData(week1 = 'GW10', location1 = 'Lumbar', v2 = 45000)
datas_cells[['GW10_T']] = scProcessData(week1 = 'GW10', location1 = 'Thoracic', v2 = 40000)
datas_cells[['GW11_C']] = scProcessData(week1 = 'GW11', location1 = 'Cervical', v2 = 50000, h2 = 0.05)
datas_cells[['GW11_L']] = scProcessData(week1 = 'GW11', location1 = 'Lumbar', v1 = 2000, v2 = 18000, h2 = 0.05, h1 = 0)
datas_cells[['GW11_T']] = scProcessData(week1 = 'GW11', location1 = 'Thoracic', v2 = 12000, h2 = 0.10, h1 = 0, v1 = 2000)
datas_cells[['GW20_C']] = scProcessData(week1 = 'GW20', location1 = 'Cervical', v2 = 30000)
datas_cells[['GW20_L']] = scProcessData(week1 = 'GW20', location1 = 'Lumbar', h2 = 0.1, v2 = 40000)
datas_cells[['GW20_T']] = scProcessData(week1 = 'GW20', location1 = 'Thoracic', h1 = 0.008, h2 = 0.08, v2 = 25000)
datas_cells[['GW23_S']] = scProcessData(week1 = 'GW23', location1 = 'Spinal_Whole', h1 = 0.008, h2 = 0.1, v2 = 28000)


datas_nucs = list()
datas_nucs[['GW8_S']] = snProcessData(week1 = 'GW8', location1 = 'Spinal_Whole')
datas_nucs[['GW10_C']] = snProcessData(week1 = 'GW10', location1 = 'Cervical', v2 = 55000, h2 = 10000)
datas_nucs[['GW10_L']] = snProcessData(week1 = 'GW10', location1 = 'Lumbar', v2 = 70000, h2 = 11000)
datas_nucs[['GW10_T']] = snProcessData(week1 = 'GW10', location1 = 'Thoracic', v2 = 55000, h2 = 10000)
datas_nucs[['GW11_C']] = snProcessData(week1 = 'GW11', location1 = 'Cervical', v2 = 15000, h2 = 5500)
datas_nucs[['GW11_L']] = snProcessData(week1 = 'GW11', location1 = 'Lumbar')
datas_nucs[['GW11_T']] = snProcessData(week1 = 'GW11', location1 = 'Thoracic', v2 = 30000)
datas_nucs[['GW20_C']] = snProcessData(week1 = 'GW20', location1 = 'Cervical', v2 = 18000)
datas_nucs[['GW20_L']] = snProcessData(week1 = 'GW20', location1 = 'Lumbar', v2 = 16000)
datas_nucs[['GW20_T']] = snProcessData(week1 = 'GW20', location1 = 'Thoracic', v2 = 12000, h2 = 5000)
datas_nucs[['GW23_S']] = snProcessData(week1 = 'GW23', location1 = 'Spinal_Whole', v2= 28000)


AMI = list()
calculateAMI = function(week1, location1){
  cat('\n Loading', week1, 'data\n')
  cat('\n Location:',location1,'\n')
  # Set directory to pull data from
  resdir1 = file.path(paste0(week1, '/', location1,'/Cells/analysis_output'))
  resdir2 = file.path(paste0(week1, '/', location1,'/Nuc/analysis_output'))
  # Load data
  seurat1 = readRDS(file.path(resdir1, '/', week1, '_filtered_with_ccAFv2.rds')

  seurat2 = readRDS(file.path(resdir2, '/', week1, '_filtered_with_ccAFv2.rds')



# AMI
subset1 = seurat2[,seurat2$ccAFv2 != 'Unknown']
predlab = subset1$ccAFv2
predlab = list(data.frame(predlab)$predlab)[[1]]
truelab = subset1$Phase
truelab = list(data.frame(truelab)$truelab)[[1]]
AMI1[[datas1]][[tag1]] <- round(AMI(truelab, predlab), 2)



###############################
# GSE155121 for spatial - GW4
##############################

setwd('/files/ccNN/testData/GSE155121')
# Convert h5ad to h5seurat
#Convert("gw4_all.h5ad", "gw4_all.h5seurat")
# Load data
data <- LoadH5Seurat("gw4_all.h5seurat")

seurat1 = data
# Perform QC
# Genes as gene symbols
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")
nCount = colSums(x = seurat1, slot = "counts")  # nCount_RNA
seurat1 <- AddMetaData(seurat1, nCount, col.name = "nCount_RNA")

v1 <- 3000
v2 <- 50000
h1 <- 0.01
h2 <- 0.12

pdf(file.path('GW4_all_QC_plot_to_choose_cutoffs.pdf'))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = 'nCount_RNA', ylab = 'nFeature_RNA', pch = 20)
abline(v = v1, col = 'red', lwd =3, lty =2)
text(v1,0,as.character(v1), cex = 0.75, pos = 1)
abline(v = v2, col = 'red', lwd =3, lty =2)
text(v2,0,as.character(v2), cex = 0.75, pos = 1)
abline(h = h1 , col = 'red', lwd =3, lty =2)
text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
abline (h = h2, col = 'red', lwd =3, lty =2)
text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
print(VlnPlot(seurat1, features = c("nCount_RNA", "percent.mito"), ncol = 2))
dev.off()

# Filter based off QC plot
keep.detect = which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 & seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)
seurat1 = subset(seurat1, cells=colnames(seurat1)[keep.detect])
cat('\n  Filtered to', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')

# Apply ccAFv2
seurat1 = PredictCellCycle(seurat1, gene_id='symbol')
write.csv(seurat1$ccAFv2, file.path('GW4_all_ccAFv2_calls.csv'))

# Apply ccSeurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat1 = CellCycleScoring(seurat1, s.features = s.genes, g2m.features = g2m.genes)

#--- ccAFv2 vs. week stage stacked barplot ---#
cf <- table(seurat1$ccAFv2, seurat1$week_stage)
totals <- colSums(cf)
data.frame(totals)
cnewdf <- rbind(cf, totals)

cnewdf <- cnewdf[apply(cnewdf, 1, function(x) !all(x==0)),]
cnewdf = data.frame(cnewdf)
cf_1 = matrix(ncol=length(unique(seurat1$week_stage)), nrow=length(unique(seurat1$ccAFv2)))
for(i in c(1:length(unique(seurat1$week_stage)))){
  for(n in c(1:length(unique(seurat1$ccAFv2)))) {
    cf_1[n,i] = cnewdf[n,i]/cnewdf[length(unique(seurat1$ccAFv2))+1, i]
  }
}
colnames(cf_1) = colnames(cf)
rownames(cf_1) = rownames(cnewdf)[1:7]
sub1 = rownames(data.frame(ccAFv2_colors)) %in% rownames(cf_1)
pdf(paste0('GW4_all_ccAFv2_percentages_020624.pdf'))
par(mar = c(8, 8, 8, 8) + 2.0)
barplot(cf_1, xlab = "", ylab = "Cell Percentage", las=2, legend.text = rownames(cf_1),  col = ccAFv2_colors[sub1], args.legend=list(x=ncol(cf_1) + 15, y=max(colSums(cf_1)), bty = "n"))
dev.off()


###############################
# U5 FUCCI G1
##############################

setwd('/files')

# Load in U5 ccAFv2 calls
U5_ccAFv2 = read.csv('ccNN/U5_ccAFv2_calls_020224.csv', row.names = 'X')
df = data.frame(table(U5_ccAFv2))
df$total = (df$Freq/dim(U5_ccAFv2)[1])*100

# Load in FUCCI U5 G1 ccAFv2 calls
U5_G1_ccAFv2 = read.csv('ccNN/testData/U5/U5_G1/analysis_output/U5_G1_ccAFv2_calls.csv', row.names = 'X')
df2 = data.frame(table(U5_G1_ccAFv2))
df2$total = (df2$Freq/dim(U5_G1_ccAFv2)[1])*100

# Calculate the fold change between ccAFv2 calls
fc = df2$total/df$total
log2fc = log2(fc)
df3 = data.frame(df$x, log2fc)
