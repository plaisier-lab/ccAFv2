##########################################################
## ccAFv2: Testing ccAFv2 models                        ##
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

#------------------------
# Set up / imports
#-----------------------
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(SeuratWrappers)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library(aricode)
use_python('/usr/bin/python3')

# Set working directory
#setwd("files/")
resdir = 'data'
savedir = 'results/ccAFv2'
dir.create(file.path(savedir, 'figures'), showWarnings = FALSE)
dir.create(file.path(savedir, 'metrics'), showWarnings = FALSE)
figdir = file.path(savedir, 'figures')
datadir = file.path(savedir, 'metrics')

# colors and order for plotting
ccSeurat_order = c('G1', 'S', 'G2M')
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")
ccAF_colors = c("G1/other" = "#9aca3c", "Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca")
ccAF_order = c("G1/other", 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1')
ccAFv2_order = c('G1/other','Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("G1/other"= "#9aca3c","Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")

# Functions
classifyCellCycle = function(seurat0, cutoff=0.5, assay='SCT', species='human', gene_id='ensembl') {
    cat('Running ccAFv2:\n')
    # Make a copy of object
    seurat1 = seurat0
    # Run SCTransform on data being sure to include the mgenes
    if(assay=='SCT') {
        seurat1 = SCTransform(seurat1, return.only.var.genes=FALSE, verbose=FALSE)
    }
    # Subset data marker genes to marker genes included in classification
    sub_genes = intersect(row.names(seurat1),mgenes)
    seurat_subset = subset(seurat1, features = sub_genes)
    # Find missing genes and assign 0s to each cell
    cat(paste0('  Total possible marker genes for this classifier: ', length(mgenes),'\n'))
    if(assay=='SCT') {
        missing_genes = setdiff(mgenes, rownames(seurat_subset[[assay]]@scale.data))
        cat(paste0('    Marker genes present in this dataset: ', nrow(seurat_subset[[assay]]@scale.data),'\n'))
        cat(paste0('    Missing marker genes in this dataset: ', length(missing_genes),'\n'))
        ## Add ERROR later to warn if not enough marker genes ##
        input_mat = seurat_subset[[assay]]@scale.data
    } else {
        missing_genes = setdiff(mgenes, rownames(seurat_subset[[assay]]@data))
        cat(paste0('    Marker genes present in this dataset: ', nrow(seurat_subset[[assay]]@data),'\n'))
        cat(paste0('    Missing marker genes in this dataset: ', length(missing_genes),'\n'))
        ## Add ERROR later to warn if not enough marker genes ##
        input_mat = seurat_subset[[assay]]@data
    }
    input_mat_scaled = t(scale(t(as.matrix(input_mat))))
    tmp = matrix(min(input_mat_scaled,na.rm=T),nrow=length(missing_genes), ncol=ncol(seurat1))
    rownames(tmp) = missing_genes
    colnames(tmp) = colnames(seurat_subset)
    input_mat_scaled_add_missing_genes = rbind(input_mat_scaled, tmp)[mgenes,]
    input_mat_scaled_add_missing_genes[!is.finite(input_mat_scaled_add_missing_genes)] = 0
    cat(paste0('  Predicting cell cycle state probabilities...\n'))
    predictions1 = predict(ccAFv2, t(input_mat_scaled_add_missing_genes))
    #colnames(predictions1) = classes
    colnames(predictions1) = c("G1", "G2/M", "Late G1", "M/Early G1", "Neural G0", "S", "S/G2")
    rownames(predictions1) = colnames(seurat1)
    df1 = data.frame(predictions1)
    cat(paste0('  Choosing cell cycle state...\n'))
    CellCycleState = data.frame(factor(colnames(predictions1)[apply(predictions1,1,which.max)], levels=c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1','Unknown')), row.names = rownames(predictions1))
    colnames(CellCycleState) = 'ccAFv2'
    df1[,'ccAFv2'] = CellCycleState$ccAFv2
    df1[which(apply(predictions1,1,max)<cutoff),'ccAFv2'] = 'Unknown'
    cat('  Adding probabilitities and predictions to metadata\n')
    seurat0 = AddMetaData(object = seurat0, metadata = df1)
    cat('Done\n')
    return(seurat0)
}

# Load ccAFv2 models and marker genes
# layers
#tags = c('700_400', '700_300', '700_200', '700_100', '600_400', '600_300', '600_200', '600_100', '500_400', '500_300', '500_200', '500_100', '400_300', '400_200', '400_100', '300_200', '300_100', '200_100')
tags = c('600_200')

# Load in data
AMI1 = list()
CellsPredict = list()
#datas = c('U5', 'LGG275_GF', 'BT322', 'BT324', 'BT326', 'BT333', 'BT363', 'BT368', 'BT363_tumor', 'BT368_tumor')
datas = c('U5', 'LGG275_GF', 'BT322', 'BT324', 'BT326', 'BT333', 'BT363', 'BT368')
for(datas1 in datas){
  cat('\n', datas1,'\n')
  AMI1[[datas1]] = list()
  CellsPredict[[datas1]] = list()
  if(datas1 == 'U5'){
    #resdir2 = file.path(resdir, datas1, 'analysis_output')
    resdir2 = file.path(resdir, datas1)
    resdir3 = file.path(resdir, datas1, 'seurat_objects')
    seurat2 = readRDS(file.path(paste0(resdir,'/', datas1,'/', datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 'X')
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if (datas1 == 'LGG275_GF'){
    #resdir2 = file.path(paste0(resdir,'/LGG/', strsplit(datas1, split = "_")[[1]][1], '/', datas1, '/analysis_output'))
    resdir2 = file.path(paste0(resdir,'/LGG/', strsplit(datas1, split = "_")[[1]][1], '/', datas1))
    resdir3 = file.path(paste0(resdir,'/LGG/', strsplit(datas1, split = "_")[[1]][1], '/', datas1, '/seurat_objects'))
    seurat2 = readRDS(file.path(resdir2, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 'X')
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
    ccAF_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccAF_calls.csv')), row.names = 'X')
    seurat2 = AddMetaData(seurat2, metadata = ccAF_calls$ccAF, col.name = 'ccAF')
    # Downstream analyses
    seurat2 <- RunPCA(seurat2, dims = 1:30, verbose=FALSE)
    seurat2 <- FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
    seurat2 <- FindClusters(seurat2, verbose=FALSE, resolution = 0.61)
    seurat2 <- RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  }else if(datas1 %in% c('BT322', 'BT324', 'BT326', 'BT333', 'BT363', 'BT368')){
    resdir2 = file.path(resdir,'GSC', datas1)
    seurat2 = readRDS(file.path(resdir2, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 1)
    colnames(phase_calls) = c('x')
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if(datas1 %in% c('BT363_tumor', 'BT368_tumor')){
    resdir = file.path(paste0('testData/GBM_tumors/',datas1, '/seurat_objects'))
    resdir2 = file.path(paste0('testData/GBM_tumors/',datas1, '/analysis_output'))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_normalized_ensembl.rds')))
    phase_calls = read.csv(file.path(resdir2, paste0(datas1, '_ccSeurat_calls.csv')), row.names = 1)
    seurat2 = AddMetaData(seurat2, metadata = phase_calls$x, col.name = 'Phase')
  }else if(datas1 == 'ScienCell_hBMMSC'){
    resdir = file.path(paste0('testData/hBMMSCs/',datas1, '/seurat_objects'))
    seurat2 = readRDS(file.path(resdir, paste0(datas1, '_processed.rds')))
  }
  for(tag1 in tags){
    cat('\n',tag1,'ccAFv2 model\n')
    # Load model and marker genes
    ccAFv2 = load_model_hdf5(file.path(paste0('model/ccAFv2_layers_', tag1,'.h5')))
    mgenes = read.csv(file.path(paste0('model/ccAFv2_genes_layers_', tag1, '.csv')))[,2]
    # Classify with ccAFv2
    seurat2 = classifyCellCycle(seurat2)
    # Calculate AMI
    subset1 = seurat2[,seurat2$ccAFv2 != 'Unknown']
    predlab = subset1$ccAFv2
    predlab = list(data.frame(predlab)$predlab)[[1]]
    truelab = subset1$Phase
    truelab = list(data.frame(truelab)$truelab)[[1]]
    AMI1[[datas1]][[tag1]] <- round(AMI(truelab, predlab), 2)
    CellsPredict[[datas1]][[tag1]] = round((dim(subset1)[2]/dim(seurat2)[2])*100, 2)
    # Plotting order & colors
    # ccAF & ccAFv2
    sub1 = ccSeurat_order %in% factor(seurat2$Phase)
    seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub1])
    sub2 = ccAF_order %in% factor(seurat2$ccAF)
    seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub2])
    sub3 = ccAFv2_order %in% factor(seurat2$ccAFv2)
    seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub3])
    p1 = DimPlot(seurat2, reduction = "umap", label=F, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub1])
    c1 = DimPlot(seurat2, reduction = "umap", label=F, pt.size = 1, group.by = 'ccAF', cols = ccAF_colors[sub2])
    c2 = DimPlot(seurat2, reduction = "umap", label=F, pt.size = 1, group.by = 'ccAFv2', cols = ccAFv2_colors[sub3])
    pdf(file.path(figdir,paste0('ccAFv2_layers_', tag1,'_', datas1, '.pdf')), width = 10, height = 8)
    # umaps
    lst = list(p1, c1, c2)
    grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,3)), top = "")
    # Stacked barplot
    # ccAFv2 vs. ccAF percentages
    df <- table(seurat2$ccAFv2, seurat2$ccAF)
    totals <- colSums(df)
    data.frame(totals)
    newdf <- rbind(df, totals)
    df_1 = matrix(ncol=length(unique(seurat2$ccAF)), nrow=length(unique(seurat2$ccAFv2)))
    for(i in c(1:length(unique(seurat2$ccAF)))){
      print(i)
      for(n in c(1:length(unique(seurat2$ccAFv2)))) {
        df_1[n,i] = newdf[n,i]/newdf[length(unique(seurat2$ccAFv2))+1, i]
      }
    }
    colnames(df_1) = colnames(df)
    rownames(df_1) = rownames(df)
    par(mar = c(8, 8, 8, 8) + 0.2)
    barplot(df_1, col = ccAFv2_colors[sub2], xlab = "", ylab = "Cell Percentage", las=2,legend.text = rownames(table(seurat2$ccAFv2, seurat2$ccAF)), args.legend=list(title="ccAFv2", x ='topright', bty='n', inset=c(-0.35,0)))
    mtext("ccAF", side=1, line=5)
    dev.off()
  }
}

# Save out AMI and Cells Predicted for each ccAFv2 model
tmp = as.data.frame(do.call(cbind, AMI1))
tmp$layers = rownames(tmp)
df = data.frame(lapply(tmp, as.character), stringsAsFactors=FALSE)
df$layers <- paste(df$layers, sep="_", 'v2')
write.csv(df, 'model/testing/FINAL/AMI_v2.csv')

tmp2 = as.data.frame(do.call(cbind, CellsPredict))
tmp2$layers = rownames(tmp2)
df2 = data.frame(lapply(tmp2, as.character), stringsAsFactors=FALSE)
df2$layers <- paste(df2$layers, sep="_", 'v2')
write.csv(df2, 'model/testing/FINAL/CellsPredicted_v2.csv')
