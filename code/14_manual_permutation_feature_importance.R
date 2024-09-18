##########################################################
## manual_permutation_feature_importance.R              ##
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


############
## Docker ##
############

# docker run -it -v '/media/old_home/home/cplaisier/ccAFv2_test:/files' cplaisier/ccafv2_seurat4


###############
## Libraries ##
###############

library(ccAFv2)
library(Seurat)
library(gtools)
library(ggplot2)


###############
## Load data ##
###############

# List of classifier genes
mgenes = read.csv(system.file("extdata", "ccAFv2_genes.csv", package = "ccAFv2"), header = TRUE, row.names = 1)[, paste0("human_ensembl")]
mgenes_sym = read.csv(system.file("extdata", "ccAFv2_genes.csv", package = "ccAFv2"), header = TRUE, row.names = 1)[, paste0("human_symbol")]

# List of cell cycle states
states1 = c("G1", "G2.M", "Late.G1", "M.Early.G1", "Neural.G0", "S", "S.G2")
states2 = c("G1", "G2/M", "Late G1", "M/Early G1", "Neural G0", "S", "S/G2")

# Load Data
setwd('/files')
seurat_obj = readRDS('U5_normalized_ensembl.rds')

# Renormalize to have all genes
seurat_obj = SCTransform(seurat_obj, return.only.var.genes=F)

# Predict cell cycle states
seurat_obj = PredictCellCycle(seurat_obj, do_sctransform=F)

# For storing data
m1 = matrix(ncol=7,nrow=861)
rownames(m1) = mgenes
colnames(m1) = states2


#################
## Permutation ##
#################

seurat_obj1 = seurat_obj
for(feature1 in mgenes) {
    print(feature1)
    # Return data to original state
    seurat_obj1@assays$SCT@scale.data = seurat_obj@assays$SCT@scale.data
    # Permute data
    seurat_obj1@assays$SCT@scale.data[feature1,] = permute(seurat_obj1@assays$SCT@scale.data[feature1,])
    # Run ccAFv2
    seurat_obj1 = PredictCellCycle(seurat_obj1, do_sctransform=F)
    # Collect differences
    for(state1 in 1:7) {
        diff1 = data.frame(diff = seurat_obj1@meta.data[,states1[state1]]-seurat_obj@meta.data[,states1[state1]],ccAF=seurat_obj1@meta.data[,'ccAF'])
        m1[feature1, states2[state1]] = as.numeric(diff1 %>% filter(ccAF==states2[state1]) %>% summarise(avg = mean(diff)))
    }
}


###################
## Write results ##
###################

m2 = m1
rownames(m2) = mgenes_sym
write.csv(m2, 'feature_importance_8_27_2024.csv')


################
## Make plots ##
################

# Make bar plots
colors1 = c(G1 = "#f37f73", `G2/M` = "#3db270", `Late G1` = "#1fb1a9", `M/Early G1` = "#6d90ca", `Neural G0` = "#d9a428", S = "#8571b2", `S/G2` = "#db7092")
plots1 = list()
pdf('feature_importance_barplots_8_27_2024.pdf')
for(state1 in 1:7) {
    df1 = data.frame(AvgDeltaLikelihood=m2[,states2[state1]], gene = rownames(m2)) %>% top_n(-20, wt=AvgDeltaLikelihood) %>% arrange(AvgDeltaLikelihood)
    p = ggplot(data=df1, aes(x=reorder(gene,-AvgDeltaLikelihood), y=AvgDeltaLikelihood)) +
        geom_bar(stat='identity', fill=colors1[states2[state1]]) +
        theme_minimal() +
        ggtitle(states2[state1]) +
        xlab('') +
        ylab('Mean Change in Likelihood')
    print(p+coord_flip())
    plots1[[states2[state1]]] = p+coord_flip()
}
dev.off()

# Make combined plot
pdf('feature_importance_barplots_combined_1col_8_27_2024.pdf', width=7.5, height=3)
grid.arrange(plots1[[states2[5]]]+theme_classic(base_size=8), plots1[[states2[1]]]+theme_classic(base_size=8), plots1[[states2[3]]]+theme_classic(base_size=8), plots1[[states2[6]]]+theme_classic(base_size=8), plots1[[states2[7]]]+theme_classic(base_size=8),plots1[[states2[2]]]+theme_classic(base_size=8),plots1[[states2[4]]]+theme_classic(base_size=8), nrow=1)
dev.off()

