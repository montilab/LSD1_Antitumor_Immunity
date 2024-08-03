#Final Filtering script
library(rvcheck)
#library(singleCellTK)
library(Seurat)
library(SingleR)
library(celldex)
library(pointillism)
library(rhdf5)
library(dplyr)

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina')

obj.seurat <- readRDS('./data/202408/20240801_SP2509.rds')

#isolate immune
Idents(obj.seurat) <- 'CellType'
DimPlot(obj.seurat)
imm <- subset(obj.seurat, idents = c("T-Lymphocyte", "B-Lymphocyte", "Myeloid", "Neutrophils")) 

#PCA and UMAP
imm <- RunPCA(imm, verbose = FALSE)
ElbowPlot(imm, ndims = 40)
imm <- FindNeighbors(imm, dims = 1:20, verbose = FALSE)
imm <- FindClusters(imm, verbose = FALSE, resolution = 0.3)
imm <- RunUMAP(imm, dims = 1:20, verbose = FALSE)

DimPlot(imm, label = TRUE, group.by = 'RNA_snn_res.0.3') + NoLegend()
DimPlot(imm, label = TRUE, group.by = 'CellType') + NoLegend()

Idents(imm) <- 'RNA_snn_res.0.3'
immmarks <- FindAllMarkers(imm, max.cells.per.ident = 200)
immmarks[immmarks$avg_log2FC >0.5 & immmarks$p_val_adj < 0.05 & immmarks$cluster == '6',]$gene

immmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 3) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(imm, features = top10$gene) + NoLegend()

#FeaturePlot
FeaturePlot(imm, features =c("KRT15","BPIFB2", "CD3D", "CD3E", "TCRG-C1", "TRDC", 
                             "GATA3", "RORA", "NKG7",
                             "CD19", "CD79A", "CD14", "FCGR3", 
                                    "C1QA", "CD74", "ITGAM", "CCR2", "KLRB1C", "CLEC4D",
                             "CD86", "ITGAX", "CD209A", "IRF8"))


Idents(imm) <- 'RNA_snn_res.0.3'
imm <- RenameIdents(imm, '0' = 'Neutrophils', '1' = 'Macrophages', '2'= 'Neutrophils',
                           '3' = 'DCs', '4' = 'T/NK', '5' = 'Monocytes',
                           '6' = 'Epithelial', '7' = 'B Cells', '8' = 'Neutrophils',
                           '9' = 'Neutrophils', '10' = 'TGD', '11' = 'ILC2')
imm$ImmIdents <- Idents(imm)
DimPlot(imm, group.by = 'ImmIdents')

imm$ImmIdents <- factor(imm$ImmIdents, levels = rev(c("T/NK", "TGD", "ILC2", "B Cells", "Neutrophils",
                                                  "Monocytes", "Macrophages", "DCs", "Epithelial")))

DotPlot(imm, features =c("KRT15","BPIFB2", "CD3D", "CD3E", "TCRG-C1", "TRDC", 
                         "GATA3", "RORA", "NKG7",
                         "CD19", "CD79A", "CD14", "FCGR3", 
                         "C1QA", "CD74", "ITGAM", "CCR2", "KLRB1C", "CLEC4D",
                         "CD86", "ITGAX", "CD209A", "IRF8")) 

saveRDS(imm, './data/202408/imm.Rds')

#DCs
Idents(imm) <- imm$ImmIdents
dcs <- subset(imm, idents = c("DCs"))

dcs <- RunPCA(dcs, verbose = FALSE)
ElbowPlot(dcs, ndims = 40)
dcs <- FindNeighbors(dcs, dims = 1:20, verbose = FALSE)
dcs <- FindClusters(dcs, verbose = FALSE, resolution = 0.3)
dcs <- RunUMAP(dcs, dims = 1:20, verbose = FALSE)
DimPlot(dcs, label = T)

Idents(dcs) <- 'RNA_snn_res.0.3'
dcmarks <- FindAllMarkers(dcs)
dcmarks[dcmarks$avg_log2FC >2 & dcmarks$p_val_adj < 0.05 & dcmarks$cluster == '3',]$gene

dcmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(dcs, features = top10$gene) + NoLegend()

FeaturePlot(dcs, features = c("CD207", #LC
                              "XCR1", #cDC1
                              "ITGAE", "CCR7", #mDC
                              "ITGAX", "ITGAM")) #cDC2
FeaturePlot(dcs, features = c("CXCL9"))

Idents(dcs) <- 'RNA_snn_res.0.3'
dcs <- RenameIdents(dcs, '0' = 'LC', '1' = 'cDC2_a', '2'= 'cDC1',
                    '3' = 'cDC2_b')
dcs$DCIdents <- Idents(dcs)
DimPlot(dcs, group.by = 'DCIdents')

VlnPlot(dcs, features = c("CXCL9"), split.by = 'treatment')

#saveRDS(dcs, './data/202408/dcs.Rds')

Idents(dcs) <- dcs$DCIdents
dcsnoLC <- subset(dcs, idents = c("LC"), invert = T)
dcsnoLC <- RunPCA(dcsnoLC, verbose = FALSE)
dcsnoLC <- RunUMAP(dcsnoLC, dims = 1:20, verbose = FALSE)
DimPlot(dcsnoLC, label = T)


#Tcells
Idents(imm) <- imm$ImmIdents
tcells <- subset(imm, idents = c("T/NK", "TGD", "ILC2"))

tcells <- RunPCA(tcells, verbose = FALSE)
ElbowPlot(tcells, ndims = 40)
tcells <- FindNeighbors(tcells, dims = 1:10, verbose = FALSE)
tcells <- FindClusters(tcells, verbose = FALSE, resolution = 0.4)
tcells <- RunUMAP(tcells, dims = 1:10, verbose = FALSE)
DimPlot(tcells, label = T)

Idents(tcells) <- 'RNA_snn_res.0.4'
tmarks <- FindAllMarkers(tcells)
tmarks[tmarks$avg_log2FC >1 & tmarks$p_val_adj < 0.05 & tmarks$cluster == '0',]$gene

tmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(tcells, features = top10$gene) + NoLegend()

FeaturePlot(tcells, features = c("CD3D", "CD3E", 
                                 "CD4", "CD8A", "CD8B",
                                 "TCRG-C1", "TRDC", "FOXP3",
                             "GATA3", "RORA", "NKG7")) 
FeaturePlot(tcells, features = c("IFNG", "CXCR3"))

Idents(tcells) <- 'RNA_snn_res.0.4'
tcells <- RenameIdents(tcells, '0' = 'Tcells', '1' = 'NK', '2'= 'ILC2',
                    '3' = 'Treg','4' = 'TGD_IL17+')
tcells$TIdents <- Idents(tcells)
DimPlot(tcells, group.by = 'TIdents')

DimPlot(tcells, split.by = 'treatment')
VlnPlot(tcells, features = c("CXCR3", "IFNG"), split.by = 'treatment')

#saveRDS(tcells, './data/202408/tcells.Rds')
