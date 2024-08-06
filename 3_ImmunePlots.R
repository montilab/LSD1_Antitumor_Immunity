library(rvcheck)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

imm <- readRDS('./data/202408/imm.Rds')
dcs <- readRDS('./data/202408/dcs.Rds')
tcells <- readRDS('./data/202408/tcells.Rds')

##
pdf('./results/202408/Immune/ImmuneUMAP.pdf', width = 5.5, height = 5)
DimPlot(imm, label = T, group.by = 'ImmIdents')
DimPlot(imm, label = T, group.by = 'RNA_snn_res.0.3') + NoLegend()
dev.off()

##
pdf('./results/202408/Immune/ImmuneMarkerDotPlot.pdf', width = 9, height = 5)
Idents(imm) <- imm$ImmIdents
dp <- DotPlot(imm, features =c("CD3D", "CD3E","CD4", "CD8A","NCR1", "TCRG-C1", "TRDC","IL17A",
                               "GATA3", "RORA", 
                               "CD79A",
                               "CLEC4D",
                               "ITGAM","CD14", "FCGR3", 
                               "C1QA", "CD74",  
                               "ITGAX","CD209A", "IRF8",
                               "KRT15","BPIFB2")) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "red")
dp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggstyle()
dev.off()

####DCs
pdf('./results/202408/Immune/DC_UMAP.pdf', width = 5.5, height = 5)
DimPlot(dcs, label = T, group.by = 'DCIdents')
DimPlot(dcs, label = T, group.by = 'RNA_snn_res.0.3') + NoLegend()
dev.off()

pdf('./results/202408/Immune/DC_UMAP_Treatment.pdf', width = 9, height = 3)
DimPlot(dcs,  group.by = 'DCIdents', split.by = 'treatment')
DimPlot(dcs, group.by = 'RNA_snn_res.0.3', split.by = 'treatment') 
dev.off()

pdf('./results/202408/Immune/DC_Heatmap.pdf', width = 5.5, height = 8)
DoHeatmap(dcs, features = top10$gene) + NoLegend()
dev.off()

#dcs_noLC
pdf('./results/202408/Immune/DCnoLC_UMAP.pdf', width = 5.5, height = 5)
DimPlot(dcsnoLC, label = T, group.by = 'DCIdents')
DimPlot(dcsnoLC, label = T, group.by = 'RNA_snn_res.0.3') + NoLegend()
dev.off()

pdf('./results/202408/Immune/DCnoLC_UMAP_Treatment.pdf', width = 7, height = 3)
DimPlot(dcsnoLC,  group.by = 'DCIdents', split.by = 'treatment')
DimPlot(dcsnoLC, group.by = 'RNA_snn_res.0.3', split.by = 'treatment') 
dev.off()

##Batf3
pdf('./results/202408/Immune/DCnoLC_BATF3UMAP_Treatment.pdf', width = 9, height = 3)
FeaturePlot(dcsnoLC, features = c("BATF3"), split.by = 'treatment')
dev.off()

dcsnoLC$batf3 <- dcsnoLC@assays$RNA@data['BATF3',]
dcsnoLC$xcr1 <- dcsnoLC@assays$RNA@data['XCR1',]
dcsnoLC$cxcl9 <- dcsnoLC@assays$RNA@data['CXCL9',]
df <- dcsnoLC@meta.data

pdf('./results/202408/Immune/BATF3_CellType.pdf', width = 6, height = 4)
ggplot2::ggplot(df, aes(x=DCIdents, y=batf3, fill=DCIdents)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="Batf3 Normalized Expression",x="CellType") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) + ggstyle()
dev.off()

pdf('./results/202408/Immune/BATF3_Treatment_BoxPlots.pdf', width = 4, height = 4)
library(ggpubr)
ggplot2::ggplot(df, aes(x=treatment, y=batf3, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="Batf3 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 1.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 2.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 1.8)  +
  scale_fill_manual(values=c("blue", "red", "green3"))
dev.off()

pdf('./results/202408/Immune/CXCL9_DC_UMAP.pdf', width = 7, height = 4)
FeaturePlot(dcsnoLC, features = c("CXCL9"), split.by = 'treatment')
dev.off()

pdf('./results/202408/Immune/CXCL9DC_Treatment_BoxPlots.pdf', width = 4, height = 4)
library(ggpubr)
ggplot2::ggplot(df, aes(x=treatment, y=cxcl9, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="Cxcl9 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3"))
dev.off()


dcs$batf3 <- dcs@assays$RNA@data['BATF3',]
dcs$xcr1 <- dcs@assays$RNA@data['XCR1',]
dcs$cxcl9 <- dcs@assays$RNA@data['CXCL9',]
dcs$cxcr3 <- dcs@assays$RNA@data['CXCR3',]
df <- dcs@meta.data
write.csv(df[,c("batf3", "xcr1", "cxcl9","cxcr3", "DCIdents", "RNA_snn_res.0.3", "treatment")], './results/202408/Immune/DC_BoxPlots.csv')


###Tcells
pdf('./results/202408/Immune/TCell_UMAP.pdf', width = 5.5, height = 5)
DimPlot(tcells, label = T, group.by = 'TIdents')
DimPlot(tcells, label = T, group.by = 'RNA_snn_res.0.4') + NoLegend()
dev.off()

pdf('./results/202408/Immune/TCell_Heatmap.pdf', width = 5.5, height = 7)
DoHeatmap(tcells, features = top10$gene) + NoLegend()
dev.off()

#
tcells$TCRBC1 <- tcells@assays$RNA@data['TRBC1',]
tcells$TCRBC2 <- tcells@assays$RNA@data['TRBC2',]
tcells$IFNG <- tcells@assays$RNA@data['IFNG',]
tcells$CXCR3 <- tcells@assays$RNA@data['CXCR3',]
df <- tcells@meta.data


pdf('./results/202408/Immune/TRBC1_Tcell_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=TCRBC1, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="TCRBC1 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

pdf('./results/202408/Immune/TRBC2_Tcell_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=TCRBC2, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="TCRBC2 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

pdf('./results/202408/Immune/CXCR3_Tcell_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=CXCR3, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="CXCR3 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

pdf('./results/202408/Immune/IFNG_Tcell_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=IFNG, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="IFNG Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

pdf('./results/202408/Immune/IFNG_CXCR3_TcellUMAP.pdf', width = 7, height = 4)
FeaturePlot(tcells, features = c("IFNG", "CXCR3"), split.by = 'treatment')
dev.off()



write.csv(df[,c("TCRBC1", "TCRBC2", "IFNG", "CXCR3", "TIdents", "treatment")], './results/202408/Immune/TCellBoxplots.csv')

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
