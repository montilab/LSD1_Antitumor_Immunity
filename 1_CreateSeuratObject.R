#Final Filtering script
library(rvcheck)
#library(singleCellTK)
library(Seurat)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(pointillism)
library(rhdf5)
#qc has already been performed. 

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina')
obj <- readRDS('./data/Kroehling_sctkQC.rds')

#first take only SP2509, control and 4nqo control
sampleswewant <- c("KM_BM_0318_1C", "KM_BM_0318_5NC", "KM_BM_0319_1", "KM_BM_0319_5", "KM_BM_0319_8", "KM_BM_0318_8NL")
obj1 <- obj[,obj$sample %in% sampleswewant]

#remove mito, decontx and doublets (intersection of all 4 algorithms)
obj2 <- obj1[ , !is.na(obj1$subsets_mito_percent) & obj1$subsets_mito_percent <= 20]
obj3 <- obj2[ , !is.na(obj2$decontX_contamination) & obj2$decontX_contamination <= .5]
a <- colnames(obj3[ ,obj3$scDblFinder_class == "doublet"])
b <- colnames(obj3[ ,obj3$scds_cxds_call == "TRUE"])
c <- colnames(obj3[ ,obj3$scds_bcds_call == "TRUE"])
d <- colnames(obj3[ ,obj3$scds_hybrid_call == "TRUE"])
com <- intersect(a,b)
com1 <- intersect(com, c)
com2 <- intersect(com1, d)
obj4 <- obj3[,! colnames(obj3) %in% com2]

#fix gene names
genes <- rownames(obj4)
genes1 <- gsub(".*_", "", genes)
genes2 <- toupper(genes1)
rownames(obj4) <- genes2


#duplicates?
n_occur <- data.frame(table(rownames(obj4)))
n_occur[n_occur$Freq > 1,]

#convert unnormalized to seurat
obj.seurat <- as.Seurat(obj4, counts = "counts", data = NULL)

#Remove other umaps made by decontx that confuse seurat
obj.seurat@reductions$decontX_KM_BM_0318_1C_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_2E_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_3E_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_4E_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_5NC_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_6E_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_7E_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0318_8NL_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_1_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_2_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_3_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_4_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_5_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_6_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_7_UMAP <- NULL
obj.seurat@reductions$decontX_KM_BM_0319_8_UMAP <- NULL

#First RNA normalize
obj.seurat <- RenameAssays(object = obj.seurat, originalexp = 'RNA')
DefaultAssay(obj.seurat) <- 'RNA'
obj.seurat <- NormalizeData(obj.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
obj.seurat <- FindVariableFeatures(obj.seurat, nfeatures = 5000, selection.method = "vst")
obj.seurat <- ScaleData(obj.seurat, features = rownames(obj.seurat))


#obj.seurat <- SCTransform(obj.seurat, assay = "RNA", method="glmGamPoi", verbose = FALSE) #vars.to.regress = "percent.mt"

#PCA and UMAP
obj.seurat <- RunPCA(obj.seurat, verbose = FALSE)
ElbowPlot(obj.seurat, ndims = 40)
obj.seurat <- FindNeighbors(obj.seurat, dims = 1:30, verbose = FALSE)
obj.seurat <- FindClusters(obj.seurat, verbose = FALSE)
obj.seurat <- RunUMAP(obj.seurat, dims = 1:30, verbose = FALSE)

DimPlot(obj.seurat, label = TRUE) + NoLegend()

#cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
#grep("PCNA", rownames(obj.seurat1), ignore.case = TRUE)
g2m.genes <- cc.genes$g2m.genes
obj.seurat <- CellCycleScoring(obj.seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(obj.seurat, group.by = 'Phase')

######use singleR to classify cells
#add singleR celltype classification
#sumexp <- pointillism::as(obj.seurat, "SummarizedExperiment")
#change rownames
#genes <- rownames(sumexp)
#genes1 <- gsub(".*-", "", genes)
#rownames(sumexp) <- genes1

genemtx <- GetAssayData(obj.seurat, assay = "RNA", slot = "data")


#ImmGen reference consists of microarray profiles of pure mouse immune cells from the project of the same name (Heng et al. 2008)
ref <- ImmGenData()
rownames(ref) <- toupper(rownames(ref))
labels <- ref$label.main
predicted <- SingleR(test = genemtx, ref = ref, assay.type.test=1,
                     labels = labels)
obj.seurat$singleRident <- predicted$labels
obj.seurat$singleRcelltype <- predicted$pruned.labels
DimPlot(obj.seurat, group.by = "singleRident", label = TRUE)
DimPlot(obj.seurat, group.by = "singleRcelltype", label = TRUE)

#rename idents
Idents(obj.seurat) <- 'RNA_snn_res.0.8'
a <- DimPlot(obj.seurat, label= T) + NoLegend()
b <- DimPlot(obj.seurat, group.by = 'singleRcelltype', label = T) + NoLegend()
a+b
obj.seurat <- RenameIdents(obj.seurat, '0' = 'Fibroblasts', '1' = 'Epithelial', '2'= 'Epithelial',
                     '3' = 'Epithelial', '4' = 'Endothelial', '5' = 'Epithelial',
                     '6' = 'Epithelial', '7' = 'Epithelial', '8' = 'Endothelial',
                     '9' = 'Fibroblasts', '10' = 'Myeloid and B', '11' = 'Endothelial',
                     '12' = 'Epithelial', '13' = 'Epithelial', '14' = 'Epithelial',
                     '15' = 'Endothelial', '16' = 'RBC', '17' = 'Endothelial', #check 16
                     '18' = 'Epithelial', '19' = 'GlialSchwann', '20' = 'Endothelial',           #check 20
                     '21' = 'Neutrophils', '22' = 'Neutrophils', '23' = 'Fibroblasts',
                     '24' = 'StromalVSMC', '25' = 'T Lymphoid and B', '26' = 'Epithelial',        #check 24
                     '27'='GlialSchwann', '28'='Lymphatic Endothelial', '29'='Epithelial',       #check 27, 28
                     '30' = 'Neutrophils', '31' = 'Fibroblasts', '32' = 'Fibroblasts',
                     '33' = 'RBC')
obj.seurat$CellType <- Idents(obj.seurat)
DimPlot(obj.seurat, group.by = 'CellType', label = T)

DimPlot(obj.seurat, group.by = 'CellType', split.by = 'treatment')

#
marks <- FindAllMarkers(obj.seurat, max.cells.per.ident = 500)
marks <- FindAllMarkers(obj.seurat, test.use = 'MAST', latent.vars = 'batch')
marks[marks$avg_log2FC > 2 & marks$p_val_adj < 0.05 & marks$cluster == 'Stromal',]$gene
##do this without max cells, with MAST and write.csv

obj.seurat$CellType <- factor(obj.seurat$CellType, levels = c("Epithelial", "Endothelial", "Lymphatic Endothelial",
                                                              "Fibroblasts",
                                                              "StromalVSMC", "GlialSchwann", "RBC", "T Lymphoid and B", "Myeloid and B",
                                                              "Neutrophils"))
obj.seurat$treatment <- facter(obj.seurat$treatment, levels = c("control", "4NQO control", "4NQO+ SP2509"))

#dotplot of markers
library(RColorBrewer)
library(ggplot2)
DotPlot(obj.seurat, features =c("KRT14", "KRT5",
                                "PECAM1","NRP1", "PDPN",
                                "COL1A1", "ACTA2", "MYH11",
                                "CD3D", "CD19", "CD14", "FCGR3", 
                                "CD74", "ITGAM",
                                "PLP1", "MBP", "HBB-BT")) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
#FeaturePlot
FeaturePlot(obj.seurat, features =c("KRT14", "KRT5",
                                "PECAM1","NRP1", "PDPN",
                                "COL1A1",
                                "CD3D", "CD19", "CD14", "FCGR3", 
                                "CD74", "ITGAM"))



saveRDS(obj.seurat, './data/202408/20240801_SP2509.rds')

