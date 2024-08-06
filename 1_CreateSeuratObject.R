#Final Filtering script
library(rvcheck)
#library(singleCellTK)
library(Seurat)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(pointillism)
library(rhdf5)
library(ggplot2)
library(RColorBrewer)
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
obj.seurat <- FindClusters(obj.seurat, verbose = FALSE, resolution = 1)
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

##add batches
obj.seurat$batches <- with(obj.seurat, ifelse((obj.seurat$sample == "KM_BM_0318_1C" | obj.seurat$sample == "KM_BM_0318_2E" | obj.seurat$sample == "KM_BM_0318_3E" | obj.seurat$sample == "KM_BM_0318_4E" | obj.seurat$sample == "KM_BM_0318_5NC" | obj.seurat$sample == "KM_BM_0318_6E" | obj.seurat$sample == "KM_BM_0318_7E" | obj.seurat$sample == "KM_BM_0318_8NL"), "batch18", "batch19"))
DimPlot(obj.seurat, split.by = "batches")

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
Idents(obj.seurat) <- 'RNA_snn_res.1'
a <- DimPlot(obj.seurat, label= T) + NoLegend()
b <- DimPlot(obj.seurat, group.by = 'singleRcelltype', label = T) + NoLegend()
a+b
obj.seurat <- RenameIdents(obj.seurat, '0' = 'Fibroblasts', '1' = 'Epithelial', '2'= 'Endothelial',
                     '3' = 'Endothelial', '4' = 'Epithelial', '5' = 'Epithelial',
                     '6' = 'Epithelial', '7' = 'Epithelial', '8' = 'Epithelial',
                     '9' = 'Fibroblasts', '10' = 'Epithelial', '11' = 'Epithelial',
                     '12' = 'Epithelial', '13' = 'Epithelial', '14' = 'Endothelial',
                     '15' = 'Myeloid', '16' = 'Endothelial', '17' = 'Epithelial', 
                     '18' = 'RBCs', '19' = 'Glial', '20' = 'Endothelial',           
                     '21' = 'Neutrophils', '22' = 'Neutrophils', '23' = 'Epithelial',
                     '24' = 'Epithelial', '25' = 'Epithelial', '26' = 'Fibroblasts',       
                     '27'='T-Lymphocyte', '28'='Glial', '29'='Lymphatic Endothelial',       
                     '30' = 'StromalVSMC', '31' = 'B-Lymphocyte', '32' = 'Epithelial',
                     '33' = 'Epithelial', '34' = 'StromalVSMC', '35' = 'Neutrophils', '36' = 'Fibroblasts', '37' = 'Fibroblasts')
#####refine stromal cluster?/ break up?
obj.seurat$CellType <- Idents(obj.seurat)
DimPlot(obj.seurat, group.by = 'CellType', label = T)

DimPlot(obj.seurat, group.by = 'CellType', split.by = 'treatment')

#
Idents(obj.seurat) <- 'RNA_snn_res.1'
#marks <- FindAllMarkers(obj.seurat, max.cells.per.ident = 500) #don't use
marks <- FindAllMarkers(obj.seurat, test.use = 'MAST', latent.vars = 'batches', only.pos = T, logfc.threshold = 0.5)
write.csv(marks, './results/202408/wholeUMAP/ClusterTypeMarkersMAST.csv')
Idents(obj.seurat) <- 'CellType'
marks <- FindAllMarkers(obj.seurat, test.use = 'MAST', latent.vars = 'batches', only.pos = T, logfc.threshold = 0.5)
write.csv(marks, './results/202408/wholeUMAP/CellTypeMarkersMAST.csv')

marks[marks$avg_log2FC > 1 & marks$p_val_adj < 0.05 & marks$cluster == '20',]$gene
##do this without max cells, with MAST and write.csv

obj.seurat$CellType <- factor(obj.seurat$CellType, levels = rev(c("Epithelial", "Endothelial", "Lymphatic Endothelial",
                                                              "Fibroblasts",
                                                              "StromalVSMC", "Glial", "RBCs",
                                                              "T-Lymphocyte", "B-Lymphocyte","Myeloid","Neutrophils")))
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
                                    "COL1A1", "ACTA2", "MYH11",
                                    "CD3D", "CD19", "CD14", "FCGR3", 
                                    "CD74", "ITGAM",
                                    "PLP1", "MBP", "HBB-BT"))


#saveRDS(obj.seurat, './data/202408/20240801_SP2509.rds')


#Figures
#

##whole umap
pdf('./results/202408/wholeUMAP/WholeUMAP.pdf', width = 6.5, height = 5)
DimPlot(obj.seurat, group.by = 'CellType', label = T, repel = T, label.size = 5) 
dev.off()

pdf('./results/202408/wholeUMAP/WholeUMAP_Cluster.pdf', width = 6.5, height = 5)
DimPlot(obj.seurat, group.by = 'RNA_snn_res.1', label = T, label.size = 5) 
dev.off()

pdf('./results/202408/wholeUMAP/DotPlot_CellTypeMarkers.pdf', width = 9, height = 5)
Idents(obj.seurat) <- 'CellType'
dp <- DotPlot(obj.seurat, features =c("KRT14", "KRT5",
                                "PECAM1","NRP1", "PDPN",
                                "COL1A1", "ACTA2", "MYH11",
                                "PLP1", "MBP", "HBB-BT",
                                "CD3D", "CD79A","CD74",  "CD14", "FCGR3")) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "red")
#+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
dp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggstyle()
dev.off()

pdf('./results/202408/wholeUMAP/FeaturePlot_CellTypeMarkers.pdf', width = 9, height = 9)
FeaturePlot(obj.seurat, features =c("KRT14", "KRT5",
                                      "PECAM1","NRP1", "PDPN",
                                      "COL1A1", "ACTA2", "MYH11",
                                      "PLP1", "MBP", "HBB-BT",
                                      "CD3D", "CD79A","CD74",  "CD14", "FCGR3"))
dev.off()


###KDM1A
pdf('./results/202408/wholeUMAP/KDM1A.pdf', width = 5, height = 4)
FeaturePlot(obj.seurat, features = c("KDM1A"))
dev.off()

pdf('./results/202408/wholeUMAP/KDM1A_bytreatment.pdf', width = 8.5, height = 3.5)
FeaturePlot(obj.seurat, features = c("KDM1A"), split.by = 'treatment')
dev.off()

obj.seurat$KDM1A <- obj.seurat@assays$RNA@data['KDM1A',]
df <- obj.seurat@meta.data

pdf('./results/202408/wholeUMAP/KDM1A_CellTypeboxplot.pdf', width = 7.5, height = 7)
ggplot2::ggplot(df, aes(x=CellType, y=KDM1A, fill=CellType)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="Normalized KDM1A Expression",x="CellType") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  ggstyle() + NoLegend()
dev.off()

pdf('./results/202408/wholeUMAP/KDM1A_CellTypeboxplot_withTreatment.pdf', width = 12, height = 12)
ggplot2::ggplot(df, aes(x=treatment, y=KDM1A, fill=treatment)) +
  geom_boxplot(width=0.7) +
  facet_wrap(facet = 'CellType') +
  labs(title="", y="Normalized KDM1A Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 2.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 3.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 2.8)  + #ylim(0,2.6) +
  scale_fill_manual(values=c("yellow", "grey", "red")) + ggstyle()
dev.off()
write.csv(df[,c("KDM1A", "RNA_snn_res.1", "CellType", "treatment")], './results/202408/wholeUMAP/KDM1A_CellTypeboxplot_withTreatment.csv')

Idents(obj.seurat) <- 'CellType'
celltypemarks <- FindAllMarkers(obj.seurat, test.use = 'MAST', latent.vars = 'batches')
celltypemarks <- FindAllMarkers(obj.seurat)
write.csv(celltypemarks, './results/202408/wholeUMAP/CellTypeMarkers.csv')

Idents(obj.seurat) <- 'RNA_snn_res.1'
clustermarks <- FindAllMarkers(obj.seurat, test.use = 'MAST', latent.vars = 'batches')
write.csv(clustermarks, './results/202408/wholeUMAP/ClusterTypeMarkers.csv')

pdf('./results/202408/wholeUMAP/H2K1_H2D1_FeaturePlot.pdf', width = 7, height = 4)
FeaturePlot(obj.seurat, features =c("H2-K1", "H2-D1"), split.by = 'treatment')
dev.off()

epi <- subset(obj.seurat, idents = c("Epithelial"))
pdf('./results/202408/Epi/H2K1_H2D1_FeaturePlot_Epi.pdf', width = 7, height = 4)
FeaturePlot(epi, features =c("H2-K1", "H2-D1"), split.by = 'treatment')
dev.off()

ggstyle <- function(font = "Helvetica", scale = 1) {
  fs <- function(x) x * scale # Dynamic font scaling
  ggplot2::theme(
    plot.title = ggplot2::element_text(family = font, size = fs(26), face = "bold", color = "#222222"),
    plot.subtitle = ggplot2::element_text(family = font, size = fs(18), margin = ggplot2::margin(0, 0, 5, 0)),
    plot.caption = ggplot2::element_blank(),
    legend.position = "right",
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.title = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.text = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
    # axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "#222222"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(size = fs(22), hjust = 0)
  )
}
