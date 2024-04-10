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

#docker run -it -v '/home/soconnor/old_home/ccNN/ccAFv2:/files' cplaisier/ccafv2_extra

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
#setwd("files/")

# Load ccSeurat phase gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
vec<-unlist(cc.genes)
cc_genes_list = data.frame(vec)$vec

# Cyclone functions and data necessary
load("compare_classifiers/pairs_functions.RData")
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

# Functions
abs_log <- function(x){
  x[x==0] <- 1
  si <- sign(x)
  si * log2(si*x)
}

# Load ccAFv2
mgenes = read.csv(system.file("extdata", "ccAFv2_genes.csv", package = "ccAFv2"), header = TRUE, row.names = 1)[, paste0('human_ensembl')]

#------------------------------------------------------
# Set up directory structure / folders
#---------------------------------------------------

tag = 'BT322'
resdir = 'data'
resdir2 = file.path(resdir, 'GSC', tag)
savedir = 'compare_classifiers'

#------------------------------------------------------
# Cross validation
#---------------------------------------------------

# Set parameters for analysis
nfolds = 1
ncores = 10

# Classifiers to compare
classifiers = c('ccafv2', 'tricycle', 'ccschwabe', 'recat', 'cyclone', 'peco')

# Create new folders for CV results
truelab = read.csv(file.path(resdir2, 'BT322_ccSeurat_calls.csv'), row.names = 'X')
for(class1 in classifiers){
  cat('\n Classifier:', toupper(class1),'\n')
  # Create new folders for CV results
  dir.create(file.path(savedir, class1), showWarnings = FALSE)
  cat('\n Loading data \n')
  if(class1 == 'ccafv2'){
    # CCAFV2
    data1 = readRDS(file.path(resdir2, paste0(tag, '_normalized_ensembl.rds')))
    data1$Phase = truelab$x
    # Save out ccAF vs. phase as csv
    #tmp1 = data.frame(data1$Phase)
    #colnames(tmp1) = 'Phase'
    #tmp2 = data.frame(data1$ccAF)
    #colnames(tmp2) = 'ccAF'
    #results = data.frame(tmp1, tmp2)
    #write.csv(results, file.path(savedir, 'ccaf/CV_classification_report_with_Phase_as_ref.csv'))
  } else if(class1 %in% c('tricycle', 'ccschwabe')){
    # TRICYCLE OR SCHWABE
    # tricycle requires gene symbols and data to be log normalized
    data1 = readRDS(file.path(resdir2, paste0(tag, '_filtered_gene_symbols.rds')))
    data1$Phase = truelab$x
    data1 = NormalizeData(data1)
    data1 = as.SingleCellExperiment(data1)
  } else if(class1 == 'peco'){
    # PECO
    # peco requires ensembl IDs and data to be quantile normalized
    data1 = readRDS(file.path(resdir2, paste0(tag, '_filtered_ensembl.rds')))
    data1$Phase = truelab$x
    # subset to 101 significant cyclic genes
    cat(' Subset to', dim(data1)[1], 'genes \n')
    data1 = subset(data1, features = rownames(training_human$sigma))
    data1 = as.SingleCellExperiment(data1)
    data1 = data_transform_quantile(data1)
  } else if(class1 == 'recat'){
    # RECAT
    data1 = readRDS(file.path(resdir2, paste0(tag, '_normalized_gene_symbols.rds')))
    data1$Phase = truelab$x
    ccAF_calls = data.frame(data1$Phase)
    data1 = GetAssayData(object = data1, slot = "counts")
    # Normalize by log2(TPM+1)
    data1 = abs_log(NormalizeTPM(data1))
    # load R scripts
    load('compare_classifiers/ola_mES_2i.RData')
    setwd('compare_classifiers/reCAT-master/R/')
    source('get_test_exp.R')
    recat_mouse_gene_symbols = mapIds(org.Mm.eg.db, keys = colnames(test_exp), keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
    recat_human_gene_symbols = convert_mouse_to_human(recat_mouse_gene_symbols)
    data1 = data1[rownames(data1) %in% recat_human_gene_symbols,]
    data1 = get_test_exp(data1)
    data1 = t(data1)
    cat(' Subset to', length(rownames(data1)), 'genes \n')
    source("get_score.R")
    setwd("../../../")
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
  write.csv(results, file.path(savedir, class1, 'BT322_CV_classification_report_with_Phase_as_ref.csv'))
}
