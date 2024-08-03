library(rvcheck)
library(Seurat)
library(pointillism)
library(rhdf5)
library(dplyr)

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina')

obj.seurat <- readRDS('./data/202408/20240801_SP2509.rds')

#isolate immune
Idents(obj.seurat) <- 'CellType'
DimPlot(obj.seurat)
epi <- subset(obj.seurat, idents = c("Epithelial")) 

#PCA and UMAP
epi <- RunPCA(epi, verbose = FALSE)
ElbowPlot(epi, ndims = 40)
epi <- FindNeighbors(epi, dims = 1:30, verbose = FALSE)
epi <- FindClusters(epi, verbose = FALSE, resolution = 0.3)
epi <- RunUMAP(epi, dims = 1:30, verbose = FALSE)

DimPlot(epi, label = TRUE, group.by = 'RNA_snn_res.0.3') + NoLegend()
DimPlot(epi, label = TRUE, group.by = 'CellType') + NoLegend()

epi$H2K1 <- epi@assays$RNA@data['H2-K1',]
epi$H2D1 <- epi@assays$RNA@data['H2-D1',]
df <- epi@meta.data

pdf('./results/202408/Epi/H2K1_Epi_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=H2K1, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="H2K1 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

pdf('./results/202408/Epi/H2D1_Epi_Treatment_BoxPlots.pdf', width = 7, height = 8)
ggplot2::ggplot(df, aes(x=treatment, y=H2D1, fill=treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="", y="H2D1 Normalized Expression",x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = list(c("4NQO control", "control")), label.y = 3.8) + 
  stat_compare_means(comparisons = list(c("4NQO+ SP2509", "control")), label.y = 4.1) +
  stat_compare_means(comparisons = list(c("4NQO control", "4NQO+ SP2509")), label.y = 3.8)  +
  scale_fill_manual(values=c("blue", "red", "green3")) + ggstyle()
dev.off()

write.csv(df[,c("H2K1", "H2D1", "treatment")], './results/202408/Epi/EpiBoxplots.csv')

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
