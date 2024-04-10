##########################################################
## ccAFv2:  CV compare cell cycle classifiers on BT322  ##
##          with ccSeurat (Phase) as reference          ##
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

#docker run -it -v '/home/soconnor/old_home/ccNN:/files' cplaisier/ccafv2_extra

#--------------------------------
# Set up section / load packages
#--------------------------------

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
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(aricode)
library(tricycle)
library(peco)
library(doParallel)
library(foreach)
library(mclust)
library(cluster)
library(TSP)
library(scran)
library(yarrr)
use_python('/usr/bin/python3')

setRepositories()
install.packages('ADImpute')
library(ADImpute)

devtools::install_github("plaisier-lab/ccafv2_R/ccAFv2")
library(ccAFv2)

# Set working directory
setwd("files/")

# Classifiers to compare
classifiers = c('ccafv2', 'tricycle', 'ccschwabe', 'recat', 'cyclone', 'peco')

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
vec<-unlist(cc.genes)
cc_genes_list = data.frame(vec)$vec

# Cyclone functions and data necessary
load("pairs_functions.RData")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

# Convert to mouse gene symbols
genes.cc_symbol = mapIds(org.Mm.eg.db, keys = genes.cc, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
genes.cc_symbol = na.omit(data.frame(genes.cc_symbol))$genes.cc_symbol

# Convert to human gene symbols
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}

genes.cc_human_symbol = convert_mouse_to_human(genes.cc_symbol)

# Load in Leng 2015 data (FUCCI sorted; a ground truth dataset)
#fucci = read.csv('GSE64016_H1andFUCCI_normalized_EC.csv', row.names = 'X')
#HumanLengESC = CreateSeuratObject(fucci)
#table(HumanLengESC$orig.ident)
#select fucci-expression cells
#colnames(HumanLengESC)[HumanLengESC$orig.ident %in% c('G1', 'S', 'G2')]
## filter them out:
#filtered_HumanLengESC <- HumanLengESC[,colnames(HumanLengESC)[HumanLengESC$orig.ident %in% c('G1', 'S', 'G2')]]

# Functions
abs_log <- function(x){
  x[x==0] <- 1
  si <- sign(x)
  si * log2(si*x)
}

# Load ccAFv2
#ccAFv2 = load_model_hdf5("model/final/ccAFv2_full_dataset_102423_fc_25_v2_layers_600_200.h5")
mgenes = read.csv("model/final/ccAFv2_genes_full_dataset_102423_fc_25_v2_layers_600_200.csv")[,2]

#------------------------------------------------------
# Set up directory structure / folders
#---------------------------------------------------

# Directory where barcodes.tsv, features.tsv, and matrix.mtx are located
resdirs = c('testData')
data_dir <- 'outs/' # where data is located
save_dir = 'analysis_output'
obj1 = 'seurat_objects'

#------------------------------------------------------
# Cross validation
#---------------------------------------------------

# Set parameters for analysis
nfolds = 10
ncores = 10

# Create new folders for CV results
analysis1 = 'cross_validation/BT322'
dir.create(analysis1, showWarnings = FALSE)
truelab = read.csv('testData/BT322_ccSeurat_phase.csv', row.names = 'X')
for(class1 in classifiers){
  cat('\n Classifier:', toupper(class1),'\n')
  dir.create(file.path(analysis1, class1), showWarnings = FALSE)
  cat('\n Loading data \n')
  if(class1 == 'ccafv2'){
    # CCAFV2
    data1 = readRDS(file.path('testData/GSC_bam/BT322/seurat_objects/BT322_normalized_ensembl.rds'))
    data1$Phase = truelab$x
  } else if(class1 %in% c('tricycle', 'ccschwabe')){
    # TRICYCLE OR SCHWABE
    # tricycle requires gene symbols and data to be log normalized
    data1 = readRDS('testData/BT322_filtered.rds')
    data1$Phase = truelab$x
    data1 = NormalizeData(data1)
    data1 = as.SingleCellExperiment(data1)
  } else if(class1 == 'peco'){
    # PECO
    # peco requires ensembl IDs and data to be quantile normalized
    # doesn't have cell cycle phase; just psuedotime - not super useful
    data1 = readRDS('testData/BT322_filtered_ensembl.rds')
    data1$Phase = truelab$x
    # subset to 101 significant cyclic genes
    cat(' Subset to', dim(data1)[1], 'genes \n')
    data1 = subset(data1, features = rownames(training_human$sigma))
    data1 = as.SingleCellExperiment(data1)
    data1 = data_transform_quantile(data1)
  } else if(class1 == 'recat'){
    # RECAT
    data1 = readRDS('testData/BT322_normalized.rds')
    data1$Phase = truelab$x
    ccAF_calls = data.frame(data1$Phase)
    data1 = GetAssayData(object = data1, slot = "counts")
    # Normalize by log2(TPM+1)
    data1 = abs_log(NormalizeTPM(data1))
    # load R scripts
    load('ola_mES_2i.RData')
    setwd('reCAT-master/R/')
    source('get_test_exp.R')
    recat_mouse_gene_symbols = mapIds(org.Mm.eg.db, keys = colnames(test_exp), keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
    recat_human_gene_symbols = convert_mouse_to_human(recat_mouse_gene_symbols)
    data1 = data1[rownames(data1) %in% recat_human_gene_symbols,]
    data1 = get_test_exp(data1)
    data1 = t(data1)
    cat(' Subset to', length(rownames(data1)), 'genes \n')
    source("get_score.R")
    setwd("../../")
  } else if(class1 == 'cyclone'){
    # CYCLONE
    data1 = readRDS('testData/BT322_filtered_ensembl.rds')
    data1$Phase = truelab$x
    data1 = NormalizeData(data1)
    indata <- which(rownames(data1) %in% unlist(hs.pairs))
    data1 <- data1[indata,] #1189 features
    cat(' Subset to', dim(data1)[1], 'genes \n')
  }
  # Initialize helper vars/indices for subsetting data (test)
  cat('\n Prepare for testing \n')
  barcodes = colnames(data1)
  nCells = dim(data1)[2]
  numSamples = round(nCells*0.8)
  allInds = 1:nCells
  sampWOreplace = replicate(nfolds,sample(allInds,numSamples,replace = FALSE))
  results = data.frame()
  testDatas = list()
  for(k in 1:nfolds){
    cat('\n Testing round', k, '\n')
    samp1 = sampWOreplace[,k]
    testDatas[[k]] = data1[,samp1]
    if(class1 == 'ccafv2'){
      # CCAFV2
      testDatas[[k]] = PredictCellCycle(testDatas[[k]])
      results = rbind(results, testDatas[[k]][[c('Phase', 'ccAFv2')]])
    } else if(class1 == 'tricycle'){
      # TRICYCLE
      testDatas[[k]] = project_cycle_space(testDatas[[k]], gname.type='SYMBOL', species='human')
      testDatas[[k]] = estimate_cycle_position(testDatas[[k]])
      testDatas[[k]] <- as.Seurat(testDatas[[k]])
      # Translate tricycle
      mylist <- list()
      for(val1 in array(testDatas[[k]]$tricyclePosition)){
        if((val1 >= 0.5*pi) & (val1 < pi)){
          mylist <- append(mylist,'S')
        } else if((val1 >= pi) & (val1 < 1.5*pi)){
          mylist <- append(mylist,'G2/M')
        } else if((val1 >= 1.5*pi) & (val1 < 1.75*pi)){
          mylist <- append(mylist,'M')
        } else if((val1 >= 1.75*pi) | (val1 < 2*pi)){
          mylist <- append(mylist,'G1/G0')
        }
      }
      # The estimated cell cycle position is bound between 0 and 2pi.
      tmp = unlist(mylist, recursive = FALSE)
      testDatas[[k]]$tricycle = tmp
      results = rbind(results, testDatas[[k]][[c('Phase', 'tricycle')]])
    } else if(class1 == 'ccschwabe'){
      # SCHWABE
      testDatas[[k]] = project_cycle_space(testDatas[[k]], gname.type='SYMBOL', species='human')
      testDatas[[k]] = estimate_Schwabe_stage(testDatas[[k]], gname.type='SYMBOL', species='human')
      testDatas[[k]] <- as.Seurat(testDatas[[k]])
      levels(testDatas[[k]]$CCStage)<-c(levels(testDatas[[k]]$CCStage),"Unknown")
      testDatas[[k]]$CCStage[is.na(testDatas[[k]]$CCStage)] <- "Unknown"
      results = rbind(results, testDatas[[k]][[c('Phase', 'CCStage')]])
    } else if(class1 == 'recat'){
      # RECAT
      score_result <- get_score(testDatas[[k]])
      recat_calls = gsub('Score', '', colnames(score_result$mean_score)[apply(score_result$mean_score, 1, which.max)])
      ccAF_calls$barcode = rownames(ccAF_calls)
      ccAF_calls_subset = ccAF_calls[rownames(ccAF_calls) %in% colnames(testDatas[[k]]),]
      ccAF_calls_subset_reordered = ccAF_calls_subset[order(factor(ccAF_calls_subset$barcode, levels = colnames(testDatas[[k]]))),]
      # colnames(testDatas[[k]]) == rownames(ccAF_calls_subset_reordered) #test
      tmp = data.frame(ccAF_calls_subset_reordered$data1.Phase, recat_calls)
      colnames(tmp) <- c('Phase','recat')
      results = rbind(results,  tmp)
    } else if(class1 == 'cyclone'){
      # CYCLONE
      results1 = cyclone(assays(as.SingleCellExperiment(testDatas[[k]]))$logcounts, pairs = hs.pairs)
      tmp = data.frame(results1$phases)
      #rownames(tmp) = colnames(testDatas[[k]])
      #colnames(tmp) = c('barcode', 'cyclone')
      tmp = data.frame(testDatas[[k]]$Phase,results1$phases)
      rownames(tmp) = colnames(testDatas[[k]])
      colnames(tmp) = c('Phase', 'cyclone')
      results = rbind(results,  tmp)
    } else if(class1 == 'peco'){
      # PECO
      pred = cycle_npreg_outsample(Y_test = testDatas[[k]], sigma_est = training_human$sigma[rownames(testDatas[[k]]),], funs_est = training_human$cellcycle_function[rownames(testDatas[[k]])], method.trend = "trendfilter", ncores=10, get_trend_estimates = FALSE)
      head(colData(pred$Y)$cellcycle_peco)
      mylist <- list()
      for(val1 in array(colData(pred$Y)$cellcycle_peco)){
        if((val1 >= 0.5*pi) & (val1 < pi)){
          mylist <- append(mylist,'S')
        } else if((val1 >= pi) & (val1 < 1.5*pi)){
          mylist <- append(mylist,'G2/M')
        } else if((val1 >= 1.5*pi) & (val1 < 1.75*pi)){
          mylist <- append(mylist,'M')
        } else if((val1 >= 1.75*pi) | (val1 < 2*pi)){
          mylist <- append(mylist,'G1/G0')
        }
      }
      # The estimated cell cycle position is bound between 0 and 2pi.
      tmp1 = unlist(mylist, recursive = FALSE)
      tmp = data.frame(testDatas[[k]]$Phase,tmp1)
      rownames(tmp) = colnames(testDatas[[k]])
      colnames(tmp) = c('Phase', 'peco')
      results = rbind(results,  tmp)
    }
  }
  cat('\n Saving out results \n')
  write.csv(results, file.path(analysis1, class1, 'CV_classification_report_020824_with_Phase_as_ref.csv'))
}


# write out ccAF calls
data1 = readRDS(file.path('testData/GSC_bam/BT322/seurat_objects/BT322_normalized_ensembl.rds'))
data1$Phase = truelab$x
tmp1 = data.frame(data1$Phase)
colnames(tmp1) = 'Phase'
tmp2 = data.frame(data1$ccAF)
colnames(tmp2) = 'ccAF'
results = data.frame(tmp1, tmp2)
write.csv(results, file.path(analysis1, 'ccaf/CV_classification_report_020824_with_Phase_as_ref.csv'))


# Apply to full dataset
# ccAFv2
data1 = PredictCellCycle(data1, gene_id = 'symbol')
write.csv(data1$ccAFv2, 'BT322_ccAFv2_calls_020824.csv')

# tricycle
# Apply to full dataset
data1 = project_cycle_space(data1, gname.type='SYMBOL', species='human')
data1 = estimate_cycle_position(data1)
data1 <- as.Seurat(data1)
# Translate tricycle
mylist <- list()
for(val1 in array(data1$tricyclePosition)){
  if((val1 >= 0.5*pi) & (val1 < pi)){
    mylist <- append(mylist,'S')
  } else if((val1 >= pi) & (val1 < 1.5*pi)){
    mylist <- append(mylist,'G2/M')
  } else if((val1 >= 1.5*pi) & (val1 < 1.75*pi)){
    mylist <- append(mylist,'M')
  } else if((val1 >= 1.75*pi) | (val1 < 2*pi)){
    mylist <- append(mylist,'G1/G0')
  }
}
tmp = unlist(mylist, recursive = FALSE)
data1$tricycle = tmp
write.csv(data1$tricycle, 'BT322_tricycle_calls_020824.csv')

# schwabe
# Apply to full dataset
data1 = project_cycle_space(data1, gname.type='SYMBOL', species='human')
data1 = estimate_Schwabe_stage(data1, gname.type='SYMBOL', species='human')
data1 <- as.Seurat(data1)
# Rename NA as unknown
levels(data1$CCStage)<-c(levels(data1$CCStage),"Unknown")
data1$CCStage[is.na(data1$CCStage)] <- "Unknown"
write.csv(data1$CCStage, 'BT322_schwabe_calls_020824.csv')

# recat
score_result <- get_score(data1)
recat_calls = gsub('Score', '', colnames(score_result$mean_score)[apply(score_result$mean_score, 1, which.max)])
df = data.frame(recat_calls)
rownames(df) = colnames(data1)
write.csv(df, 'BT322_recat_calls_020824.csv')

# cyclone
results = cyclone(assays(as.SingleCellExperiment(data1))$logcounts, pairs = hs.pairs)
tmp = data.frame(results$phases)
rownames(tmp) = colnames(data1)
write.csv(tmp, 'BT322_cyclone_calls_020824.csv')

# peco
pred <- cycle_npreg_outsample(Y_test = data1, sigma_est = training_human$sigma[rownames(data1),], funs_est = training_human$cellcycle_function[rownames(data1)], method.trend = "trendfilter",get_trend_estimates = FALSE, ncores = 10)
head(colData(pred$Y)$cellcycle_peco)
mylist <- list()
for(val1 in array(colData(pred$Y)$cellcycle_peco)){
  if((val1 >= 0.5*pi) & (val1 < pi)){
    mylist <- append(mylist,'S')
  } else if((val1 >= pi) & (val1 < 1.5*pi)){
    mylist <- append(mylist,'G2/M')
  } else if((val1 >= 1.5*pi) & (val1 < 1.75*pi)){
    mylist <- append(mylist,'M')
  } else if((val1 >= 1.75*pi) | (val1 < 2*pi)){
    mylist <- append(mylist,'G1/G0')
  }
}
# The estimated cell cycle position is bound between 0 and 2pi.
tmp = unlist(mylist, recursive = FALSE)
data1$peco = tmp
write.csv(data1$peco, 'BT322_peco_calls_020824.csv')






# Load in data
seurat2 = readRDS('testData/GSC_bam/BT322/seurat_objects/BT322_normalized_ensembl.rds')

# Load in cell cycle calls from each classifier
# ccSeurat calls
ccSeurat_calls = read.csv('BT322_ccSeurat_phase.csv', row.names = 1)
seurat2$Phase = ccSeurat_calls$x
# ccAFv2 calls
ccAFv2_calls = read.csv('BT322_ccAFv2_calls_020824.csv', row.names = 1)
seurat2$ccAFv2 = ccAFv2_calls$x
# tricycle calls
tricycle_calls = read.csv('BT322_tricycle_calls_020824.csv', row.names = 1)
seurat2$tricycle = tricycle_calls$x
# schwabe calls
schwabe_calls = read.csv('BT322_schwabe_calls_020824.csv', row.names = 1)
seurat2$schwabe = schwabe_calls$x
# recat calls
recat_calls = read.csv('BT322_recat_calls_020824.csv', row.names = 1)
seurat2$recat = recat_calls$recat_calls
# cyclone calls
cyclone_calls = read.csv('BT322_cyclone_calls_020824.csv', row.names = 1)
seurat2$cyclone = cyclone_calls$results.phases
# peco calls
peco_calls = read.csv('BT322_peco_calls_020824.csv', row.names = 1)
seurat2$peco = peco_calls$x


# Plotting order & colors
# ccAF & ccAFv2
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
sub1 = ccAFv2_order %in% factor(seurat2$ccAF)
seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub1])
sub2 = ccAFv2_order %in% factor(seurat2$ccAFv2)
seurat2$ccAFv2 <- factor(seurat2$ccAFv2, levels = ccAFv2_order[sub2])
# Seurat
ccSeurat_order = c('G1', 'S', 'G2M')
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")
sub3 = ccSeurat_order %in% factor(seurat2$Phase)
seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub3])
# Tricycle
tricycle_order = c('G1/G0', 'S', 'G2/M', 'M')
tricycle_colors = c("G1/G0" = "#f37f73", "S" = "#8571b2", "G2/M" = "#3db270", "M" = "#6d90ca", "Unknown" = "#D3D3D3")
sub4 = tricycle_order %in% factor(seurat2$tricycle)
seurat2$tricycle <- factor(seurat2$tricycle, levels = tricycle_order[sub4])
# ccSchwabe
ccSchwabe_order = c('G1.S', 'S', 'G2', 'G2.M', 'M.G1', 'Unknown')
ccSchwabe_colors = c("G1.S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2.M" = "#3db270", "M.G1"= "#6d90ca", "Unknown" = "#D3D3D3")
sub5 = ccSchwabe_order %in% factor(seurat2$schwabe)
seurat2$ccSchwabe <- factor(seurat2$schwabe, levels = ccSchwabe_order[sub5])
# reCAT
recat_order = c('G1', 'G1S', 'S', 'G2', 'G2M', 'M')
recat_colors = c("G1" = "#f37f73", "G1S" = "#1fb1a9", "S" = "#8571b2", "G2" = "#db7092", "G2M" = "#3db270", "M"= "#6d90ca")
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

# Visualize umap; adjust above parameters if need be
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
tag = 'BT322'
pdf(file.path(paste0(tag, "_visualization_test_020824.pdf")), height = 12, width = 24)
lst1 = list(c1, c2, p1, t1, r1, s1, pe1, cy1)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2, 3, 4), c(5,6,7, 8)), top = "")
dev.off()
