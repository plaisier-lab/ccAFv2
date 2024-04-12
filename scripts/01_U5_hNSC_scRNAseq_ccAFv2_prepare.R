##########################################################
## ccAFv2: U5-hNSC scRNA-seq analysis                   ##
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
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)

# Plotting order & colors
ccAFv2_order = c('Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown')
ccAFv2_colors = c("Neural G0" = "#d9a428", "G1" = "#f37f73", "Late G1" = "#1fb1a9",  "S" = "#8571b2", "S/G2" = "#db7092", "G2/M" = "#3db270" ,"M/Early G1" = "#6d90ca",  "Unknown" = "#D3D3D3")
ccSeurat_order = c('G1', 'S', 'G2M')
ccSeurat_colors = c("G1" = "#f37f73", "S" = "#8571b2", "G2M" = "#3db270")

# Gene sets
mesenchymal_geneset = c("CD44", "SPARCL1", "MKI67", "AURKA")
proneural_geneset = c("OLIG1", "OLIG2", "SOX4", "SOX8", "SOX2")
hypoxia_geneset = c("VEGFA", "ANGPTL4", "NDRG1", "CD44")
G0_G1_p53_targets = c("MDM2", "BBC3", "ZMAT3", "TP53I3", "CDKN1A")
goi_lst2 = list(mesenchymal_geneset, proneural_geneset, hypoxia_geneset, G0_G1_p53_targets)
sapply(goi_lst2, length)

neuralGO_geneset = c("SCRG1", "PLP1", "S100B", "GPM6B", "BEX1", "PTPRZ1", "PRCP", "PTN", "SOX4", "SAT1")
G1_geneset = c("IGFBP3", "IGFBP5", "MIAT", "MAP3K7CL", "AHNAK2", "TPST2", "DLG1" , "CMTM7", "C6orf15", "GJB2")
late_G1_geneset = c("EDN1", "CYR61", "ANKRD1", "CTGF", "PLK2", "UGCG", "ARID5B", "PLAU", "CCL2")
S_geneset = c("CCNE2", "CLSPN", "GINS2", "PCNA", "ATAD2", "MCM7", "MCM3", "SLBP", "GMNN", "KIAA0101")
s_g2_geneset = c("HIST1H4C", "CDK1", "HIST1H1E", "HIST1H1B", "UBE2C", "RRM2", "ZWINT", "HIST1H1C", "HMGB2")
G2_M_geneset = c("CCNB1", "CENPF", "CKS2", "PTTG1", "CDC20", "TOP2A", "NUSAP1", "CENPA")
M_early_G1_geneset = c("HMGN2", "TUBA1B", "STMN1", "BIRC5", "HMGB1", "TROAP", "HNRNPA2B1", "H2AFZ", "ARL6IP1")
goi_lst = list(neuralGO_geneset,G1_geneset, late_G1_geneset, S_geneset, s_g2_geneset, G2_M_geneset, M_early_G1_geneset )
sapply(goi_lst, length)

#------------------------------------------------------
# Read in data section / set up seurat object / QC
#---------------------------------------------------

# Set working directory
#setwd("files/")

#------------------------------------------------------
# Load in data (10X)
#---------------------------------------------------

# Set data directory as resdir
resdir = 'data'
tag = 'U5'

# Directory where barcodes.tsv, features.tsv, and matrix.mtx are located
data_dir = 'outs/filtered_feature_bc_matrix/' # where data is located

# Create folder to save out figures / csv files / etc.
save_dir = 'analysis_output'
dir.create(file.path(resdir, tag, save_dir), showWarnings = FALSE)
resdir2 = file.path(resdir, tag, save_dir)

# Create folder to save out seurat objects (unfiltered, filtered, normalized, etc.)
obj1 = 'seurat_objects'
dir.create(file.path(resdir, tag, obj1), showWarnings = FALSE)
resdir3 = file.path(resdir, tag, obj1)

# Load in data (Ensembl IDs)
data = Read10X(file.path(resdir, tag, data_dir), gene.column=1) # column 1 is ensembl IDs (in 10X mtx file)

# Create seurat object
seurat1 = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
dim(seurat1) #[1] 14734  3049

# Find mitochondrial genes
mitogenes = read_csv(file.path(resdir,'mito_genes.csv'))
seurat1[['percent.mito']] = PercentageFeatureSet(seurat1, features = as.list(mitogenes)$mito)/100

# Read in previous ccAF cell labels
ccAF = read_csv(file.path(resdir, tag, paste0(tag, '_ccAF_calls.csv')))
ccAF[ccAF == "G1/Differentiation"] = "Neural G0"
seurat1 = AddMetaData(object = seurat1, metadata = ccAF$ccAF, col.name = 'ccAF')
# Remove G1/other cells
seurat1 = subset(seurat1, subset = ccAF != 'G1/other')
dim(seurat1) #[1] 14734  2973


#----------------------
# Quality control
#----------------------

# ADJUST VERTICAL AND HORIZONTAL CUTOFF LINES DEPENDING ON DATA
v1 <- 1000
v2 <- 20000
h1 <- 0.0001
h2 <- 0.07

# Plot
pdf(file.path(resdir2, paste0(tag, "_QC_plot.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = v1, col = "red", lwd =3, lty =2)
text(v1,0,as.character(v1), cex = 0.75, pos = 1)
abline(v = v2, col = "red", lwd =3, lty =2)
text(v2,0,as.character(v2), cex = 0.75, pos = 1)
abline(h = h1 , col = "red", lwd =3, lty =2)
text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
abline (h = h2, col = "red", lwd =3, lty =2)
text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
dev.off()

# Quality control filtering
keep.detect <- which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 &
                       seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)

# Subset data
seurat1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect])
dim(seurat1) #[1] 14734  2962

# Save out filtered object as rds, h5ad, and loom
# Ensembl IDs
saveRDS(seurat1, file.path(resdir3, paste0(tag, "_filtered_ensembl.rds")))
#SaveH5Seurat(seurat1, file.path(resdir3, paste0(tag, "_filtered_ensembl.h5Seurat")), overwrite = TRUE)
#Convert(file.path(resdir3, paste0(tag, "_filtered_ensembl.h5Seurat")), dest = "h5ad", overwrite = TRUE)
#data_loom <- as.loom(seurat1, file.path(resdir3, paste0(tag, "_filtered_ensembl.loom")), verbose = FALSE, overwrite = TRUE)
#data_loom$close_all()


#---------------------------------------
# Normalization / Downstream analysis
#---------------------------------------

#seurat1 = readRDS(file.path(resdir3, paste0(tag, "_filtered_ensembl.rds")))

# Save as new object so can go back to previous non-normalized seurat obj if need to
seurat2 <- seurat1

# Normalize with SCTransform
seurat2 <- SCTransform(seurat2, verbose = FALSE)
#seurat2 <- SCTransform(seurat2, verbose = FALSE, return.only.var.genes = FALSE)
#saveRDS(seurat2, file.path(resdir3, paste0(tag, "_normalized_ensembl_all_genes.rds")))
#SaveH5Seurat(seurat2, overwrite = TRUE, file.path(resdir3, paste0(tag, "_normalized_ensembl.h5Seurat")))
#Convert(file.path(resdir3, paste0(tag, "_normalized_ensembl.h5Seurat")), overwrite = TRUE, dest = "h5ad")


# Downstream analysis
seurat2 <- RunPCA(seurat2, verbose = FALSE)
seurat2 <- FindNeighbors(seurat2, verbose = FALSE)
seurat2 <- RunUMAP(seurat2, dims=1:10, verbose = FALSE)

# Save out normalized object as rds, h5ad, and loom
saveRDS(seurat2, file.path(resdir3, paste0(tag, "_normalized_ensembl.rds")))
#SaveH5Seurat(seurat2, overwrite = TRUE, file.path(resdir3, paste0(tag, "_normalized_ensembl.h5Seurat")))
#Convert(file.path(resdir3, paste0(tag, "_normalized_ensembl.h5Seurat")), overwrite = TRUE, dest = "h5ad")
#data_loom <- as.loom(seurat3, file.path(resdir3, paste0(tag, "_normalized_ensembl.loom")), verbose = FALSE, overwrite = TRUE)
#data_loom$close_all()


#---------------------------------------
# Redo with gene symbols
#---------------------------------------
# Load in data
# Genes as gene symbols
# change genes.tsv to features.tsv
data <- Read10X(file.path(resdir, tag, data_dir), gene.column=2) # column 2 is gene symbols (in 10X mtx file)

# Substitute underscores
rownames(data) = gsub("_", "-", rownames(data))

# Create seurat object
seurat1 <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

# Read in ccAF labels from before
ccAF <- read_csv(file.path(resdir, tag, paste0(tag, "_ccAF_calls.csv")))
ccAF[ccAF == "G1/Differentiation"] <- "Neural G0"
seurat1 <- AddMetaData(object = seurat1, metadata = ccAF$ccAF, col.name = 'ccAF')
# Remove G1/other cells
seurat1 = subset(seurat1, subset = ccAF != 'G1/other')

#----------------------
# Quality control
#----------------------

#------------------------------------------------------
# Quality control
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")

# ADJUST VERTICAL AND HORIZONTAL CUTOFF LINES DEPENDING ON DATA
v1 <- 1000
v2 <- 20000
h1 <- 0.0001
h2 <- 0.07

# Quality control filtering
keep.detect <- which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 &
                       seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)

# Subset data
seurat1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect])

#----------------------------------
# Save out filtered data objects
#----------------------------------

# Save out filtered object as rds, h5ad, and loom
saveRDS(seurat1, file.path(resdir3, paste0(tag, "_filtered_gene_symbols.rds")))
#SaveH5Seurat(seurat1, file.path(resdir3, paste0(tag, "_filtered_gene_symbols.h5Seurat")), overwrite = TRUE)
#Convert(file.path(resdir3, paste0(tag, "_filtered_gene_symbols.h5Seurat")), dest = "h5ad", overwrite = TRUE)
#data_loom <- as.loom(seurat1, file.path(resdir3, paste0(tag, "_filtered_gene_symbols.loom")), verbose = FALSE, overwrite = TRUE)
#data_loom$close_all()


#------------------------------------------------------
# Normalization / Downstream analysis
#---------------------------------------------------

# Save as new object so can go back to previous non-normalized seurat obj if need to
seurat2 <- seurat1

# Normalize with SCTransform
seurat2 <- SCTransform(seurat2, verbose = FALSE)

# Classify cell cycle with CellCycleScoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat2 <- CellCycleScoring(seurat2, s.features = s.genes, g2m.features = g2m.genes)
# Save out ccSeurat classifications
write.csv(seurat2$S.Score, file.path(resdir2, paste0(tag,'_ccSeurat_S_Score.csv')))
write.csv(seurat2$G2M.Score, file.path(resdir2, paste0(tag,'_ccSeurat_G2M_Score.csv')))
write.csv(seurat2$Phase, file.path(resdir2, paste0(tag,'_ccSeurat_calls.csv')))

# Downstream analysis
seurat2 <- RunPCA(seurat2, verbose = FALSE)
seurat2 <- FindNeighbors(seurat2, verbose = FALSE)
seurat2 <- RunUMAP(seurat2, dims=1:30, verbose = FALSE)

# Save out normalized object as rds, h5ad, and loom
saveRDS(seurat2, file.path(resdir3, paste0(tag, "_normalized_gene_symbols.rds")))
#SaveH5Seurat(seurat2, overwrite = TRUE, file.path(resdir3, paste0(tag, "_normalized_gene_symbols.h5Seurat")))
#Convert(file.path(resdir3, paste0(tag, "_normalized_gene_symbols.h5Seurat")), overwrite = TRUE, dest = "h5ad")
#data_loom <- as.loom(seurat2, file.path(resdir3, paste0(tag, "_normalized_gene_symbols.loom")), verbose = FALSE, overwrite = TRUE)
#data_loom$close_all()


# Find cluster marker genes and plot
# Order cell cycle calls
sub1 = ccAFv2_order %in% factor(seurat2$ccAF)
seurat2$ccAF <- factor(seurat2$ccAF, levels = ccAFv2_order[sub1])
sub2 = ccSeurat_order %in% factor(seurat2$Phase)
seurat2$Phase <- factor(seurat2$Phase, levels = ccSeurat_order[sub2])

# Set identities to cell labels
Idents(seurat2) = seurat2$ccAF

# Find cluster marker genes
cluster_markers = FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(cluster_markers, file.path(resdir2, paste0(tag,"_scTransform_Markers_together.csv")))
cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

# Visualize UMAP
c1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6, pt.size = 1, group.by = 'ccAF', cols = ccAFv2_colors[sub1])
p1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6, pt.size = 1, group.by = 'Phase', cols = ccSeurat_colors[sub2])

# Make Violin plots for each cluster
o1 = VlnPlot(object = seurat2, features= "nFeature_RNA",   pt.size=0.1 )
o2 = VlnPlot(object = seurat2, features= "nCount_RNA",  pt.size=0.1 )
o4 = VlnPlot(object = seurat2, features= "percent.mito",   pt.size=0.1)

# Dump out pdf of plots
pdf(file.path(resdir2, paste0(tag, ".pdf")), width = 8, height = 10)
# umaps
lst = list(c1, p1)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,NA)), top = "")
lst2 = list( o1, o2, o4)
grid.arrange(grobs = lst2, layout_matrix = rbind(c(1,2), c(3)), top = "")
# Heatmap
DoHeatmap(object = seurat2, features = top10$gene)
# ccAF Heatmap
DoHeatmap(object = seurat2, features=unlist(goi_lst), group.colors = list(data.frame(ccAFv2_colors[sub1])$ccAFv2_colors.sub1.)[[1]]) + ggtitle("ccAF genes")
# Cylin ridgeplots
features <- c("CCND1", "CCNE2", "CCNA2", "CCNB1","CDK1","CDK2")
print(RidgePlot(seurat2, features, ncol=2))
# Specific gene sets = proneural, mesenchymal, hypoxia, p53 targets
p1 = lapply(goi_lst2, function(goi){
  goi = intersect(goi, rownames(seurat2))
  print(FeaturePlot(object = seurat2, features = goi, coord.fixed = TRUE))
})
dev.off()
