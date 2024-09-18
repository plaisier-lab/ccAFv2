# docker run -it -v '/media/old_home/home/cplaisier/spatial_visium_HD:/files' cplaisier/ccafv2_seurat5

# packages required for Visium HD
# install.packages("hdf5r")
# Sys.setenv(LIBARROW_MINIMAL = "false")
# Sys.setenv(ARROW_WITH_ZSTD = "ON")
# install.packages("arrow")
# remotes::install_github('satijalab/seurat-data')
# remotes::install_github('satijalab/seurat', ref='visium-hd')

# docker run -it -v '/media/old_home/home/cplaisier/spatial_visium_HD:/files' cplaisier/seurat_visium_hd

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ccAFv2)
library(reticulate)
use_python('/usr/bin/python3')

setwd('/files')

PredictCellCycle = function(seurat_obj, cutoff=0.5, do_sctransform=TRUE, assay='SCT', species='human', gene_id='ensembl', spatial = FALSE) {
    cat('Running ccAFv2:\n')
    # Make a copy of object
    seurat1 = seurat_obj

    # Load model and marker genes
    ccAFv2 = keras::load_model_hdf5(system.file('extdata', 'ccAFv2_model.h5', package='ccAFv2'))
    classes = read.csv(system.file('extdata', 'ccAFv2_classes.txt', package='ccAFv2'), header=FALSE)$V1
    mgenes = read.csv(system.file('extdata', 'ccAFv2_genes.csv', package='ccAFv2'), header=TRUE, row.names=1)[,paste0(species,'_',gene_id)]
    
    # Run SCTransform on data being sure to include the mgenes
    if(assay=='SCT' & do_sctransform) {
        cat('  Redoing SCTransform to ensure maximum overlap with classifier genes...\n')
        if(!spatial) {
            seurat1 = SCTransform(seurat1, return.only.var.genes=FALSE, verbose=FALSE)
        } else {
            seurat1 = SCTransform(seurat1, assay = 'Spatial', return.only.var.genes=FALSE, verbose=FALSE)
        }
    }

    # Subset data marker genes to marker genes included in classification
    sub_genes = intersect(row.names(seurat1),mgenes)

    # Find missing genes and assign 0s to each cell
    cat(paste0('  Total possible marker genes for this classifier: ', length(mgenes),'\n'))
    if(assay=='SCT') {
        input_mat = seurat1@assays$SCT@scale.data[sub_genes,]
    } else {
        input_mat = (seurat1@assays[[assay]])$data[sub_genes,]
    }
    missing_genes = setdiff(mgenes, rownames(input_mat))
    cat(paste0('    Marker genes present in this dataset: ', nrow(input_mat),'\n'))
    cat(paste0('    Missing marker genes in this dataset: ', length(missing_genes),'\n'))
    if(nrow(input_mat)<=689) {
        warning("Overlap below 80%: try setting 'do_sctransform' parameter to TRUE.")
    }
    input_mat_scaled = t(scale(t(as.matrix(input_mat))))
    tmp = matrix(min(input_mat_scaled,na.rm=T),nrow=length(missing_genes), ncol=ncol(seurat1))
    rownames(tmp) = missing_genes
    colnames(tmp) = colnames(input_mat)
    input_mat_scaled_add_missing_genes = rbind(input_mat_scaled, tmp)[mgenes,]
    input_mat_scaled_add_missing_genes[!is.finite(input_mat_scaled_add_missing_genes)] = 0
    cat(paste0('  Predicting cell cycle state probabilities...\n'))
    predictions1 = predict(ccAFv2, t(input_mat_scaled_add_missing_genes))
    colnames(predictions1) = classes
    rownames(predictions1) = colnames(seurat1)
    df1 = data.frame(predictions1)
    cat(paste0('  Choosing cell cycle state...\n'))
    CellCycleState = data.frame(factor(colnames(predictions1)[apply(predictions1,1,which.max)], levels=c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1','Unknown')), row.names = rownames(predictions1))
    colnames(CellCycleState) = 'ccAFv2'
    df1[,'ccAFv2'] = CellCycleState$ccAFv2
    df1[which(apply(predictions1,1,max)<cutoff),'ccAFv2'] = 'Unknown'
    cat('  Adding probabilities and predictions to metadata\n')
    seurat_obj = AddMetaData(object = seurat_obj, metadata = df1)
    cat('Done\n')
    return(seurat_obj)
}


########################
### Mouse embryo 8uM  ##
########################

mm_embryo = Load10X_Spatial(data.dir='mouse_embryo_10X', bin.size=c(8))

# Setting default assay changes between 8um and 16um binning
Assays(mm_embryo)
DefaultAssay(mm_embryo) = "Spatial.008um"

# QC
#pdf('prelim_plot.pdf')
#vln.plot = VlnPlot(mm_embryo, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
#count.plot = SpatialFeaturePlot(mm_embryo, features = "nCount_Spatial.008um") + theme(legend.position = "right")
#vln.plot | count.plot
#dev.off()

# Normalization
# mm_embryo = SCTransform(mm_embryo, assay = "Spatial.008um", verbose = FALSE)
mm_embryo = NormalizeData(mm_embryo)

# ccAFv2 spatial style
mm_embryo = PredictCellCycle(mm_embryo, species='mouse', gene_id='symbol', spatial=TRUE, do_sctransform=FALSE, assay='Spatial.008um', cutoff=0.7)

# Plot spatial style with ccAFv2
#pdf('ccAFv2_SpatialDimPlot_slice1_008um_embryo.pdf')
tiff('ccAFv2.tiff', units='in', width=5, height=5, res=1600)
SpatialDimPlot.ccAFv2(mm_embryo) + theme(legend.position = "right")
dev.off()
for(i in c('Neural.G0','G1','Late.G1','S','S.G2','G2.M','M.Early.G1')) {
    tiff(paste0(i,'.tiff'), units='in', width=5, height=5, res=1600)
    sp1 = SpatialFeaturePlot(mm_embryo, features=i) & scale_fill_gradient(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "white", high = "red4")
    print(sp1)
    dev.off()
}

for(i in c('Mki67', 'Dcx', 'Pax6','Eomes', 'Tbr1', 'Reln', 'Mef2c', 'Sstr2', 'Nr2e1', 'Wnt7b', 'Hepacam', 'Ttr', 'Krt5', 'Col1a1')) {
    tiff(paste0(i,'.tiff'), units='in', width=5, height=5, res=1600)
    sp1 = SpatialFeaturePlot(mm_embryo, features=i)
    print(sp1)
    dev.off()
}

