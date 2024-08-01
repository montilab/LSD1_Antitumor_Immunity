#Final Filtering script
#cleaned and added min nFeature cutoff 1/31/22
library(rvcheck)
library(singleCellTK)
library(Seurat)
library(glmGamPoi)
library(ggVennDiagram)
library(SingleR)
library(celldex)
library(pointillism)
library(rhdf5)
library(dplyr)
library(ggplot2)
#qc has already been performed. 

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina')
obj <- readRDS('./data/Kroehling_sctkQC.rds')
#remove mito, decontx and doublets (intersection of all 4 algorithms)
obj1 <- obj[ , !is.na(obj$subsets_mito_percent) & obj$subsets_mito_percent <= 20]
obj2 <- obj1[ , !is.na(obj1$decontX_contamination) & obj1$decontX_contamination <= .5]
a <- colnames(obj2[ ,obj2$scDblFinder_class == "doublet"])
b <- colnames(obj2[ ,obj2$scds_cxds_call == "TRUE"])
c <- colnames(obj2[ ,obj2$scds_bcds_call == "TRUE"])
d <- colnames(obj2[ ,obj2$scds_hybrid_call == "TRUE"])
com <- intersect(a,b)
com1 <- intersect(com, c)
com2 <- intersect(com1, d)
obj3 <- obj2[,! colnames(obj2) %in% com2]
obj3 <- obj3[ , obj3$detected >= 400]


genes <- rownames(obj3)
genes1 <- gsub(".*_", "", genes)
genes2 <- toupper(genes1)
rownames(obj3) <- genes2

#convert unnormalized to seurat
#no cutoff on min genes per cell, have already done that
obj.seurat <- as.Seurat(obj3, counts = "counts", data = NULL)

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
obj.seurat <- FindVariableFeatures(obj.seurat, nfeatures = 50000)
obj.seurat <- ScaleData(obj.seurat)

#sctransform calculates a model of technical noise in scRNA-seq data using 
#‘regularized negative binomial regression’. The residuals for this model 
#are normalized values, and can be positive or negative. Positive residuals 
#for a given gene in a given cell indicate that we observed more UMIs than 
#expected given the gene’s average expression in the population and cellular 
#sequencing depth, while negative residuals indicate the converse.
obj.seurat <- SCTransform(obj.seurat, assay = "RNA", method="glmGamPoi", verbose = FALSE) #vars.to.regress = "percent.mt"
obj.seurat <- RunPCA(obj.seurat, verbose = FALSE)
obj.seurat <- RunUMAP(obj.seurat, dims = 1:30, verbose = FALSE)
obj.seurat <- FindNeighbors(obj.seurat, dims = 1:30, verbose = FALSE)
obj.seurat <- FindClusters(obj.seurat, verbose = FALSE)
DimPlot(obj.seurat, label = TRUE) + NoLegend()

FeaturePlot(obj.seurat, features = "subsets_mito_percent")
DimPlot(obj.seurat, split.by = "Sample_ID", ncol = 8)

obj.seurat$batches <- with(obj.seurat, ifelse((obj.seurat$sample == "KM_BM_0318_1C" | obj.seurat$sample == "KM_BM_0318_2E" | obj.seurat$sample == "KM_BM_0318_3E" | obj.seurat$sample == "KM_BM_0318_4E" | obj.seurat$sample == "KM_BM_0318_5NC" | obj.seurat$sample == "KM_BM_0318_6E" | obj.seurat$sample == "KM_BM_0318_7E" | obj.seurat$sample == "KM_BM_0318_8NL"), "batch18", "batch19"))
DimPlot(obj.seurat, split.by = "batches")

#cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
#grep("PCNA", rownames(obj.seurat1), ignore.case = TRUE)
g2m.genes <- cc.genes$g2m.genes
obj.seurat <- CellCycleScoring(obj.seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(obj.seurat, group.by = 'Phase')

#add singleR celltype classification
sumexp <- as(obj.seurat, "SummarizedExperiment")
#change rownames
genes <- rownames(sumexp)
genes1 <- gsub(".*-", "", genes)
rownames(sumexp) <- genes1

#This reference consists of a collection of mouse bulk RNA-seq data sets downloaded from the gene expression omnibus (Benayoun et al. 2019). A variety of cell types are available, again mostly from blood but also covering several other tissues.
#ref <- MouseRNAseqData()

#ImmGen reference consists of microarray profiles of pure mouse immune cells from the project of the same name (Heng et al. 2008)
ref <- ImmGenData()
rownames(ref) <- toupper(rownames(ref))
labels <- ref$label.main
predicted <- SingleR(test = sumexp, ref = ref, assay.type.test=1,
                     labels = labels)
obj.seurat$singleRident <- predicted$labels
obj.seurat$singleRcelltype <- predicted$pruned.labels
DimPlot(obj.seurat, group.by = "singleRident", label = TRUE)
DimPlot(obj.seurat, group.by = "singleRcelltype", label = TRUE)

obj.seurat$singleRscoreEndoconf <- predicted$scores[,5]
obj.seurat$singleRscoreEpithelialconf <- predicted$scores[,7]
obj.seurat$singleRscoreFibconf <- predicted$scores[,8]
obj.seurat$singleRscoreStromalconf <- predicted$scores[,18]

FeaturePlot(obj.seurat, features = 'singleRscoreStromalconf')

#saveRDS(obj.seurat, './data/seuratobj_finalfiltering_wnFeat_01312022_a.rds')

#best object as of 02_18_2022
obj.seurat <- readRDS('./data/4nqo_filtered_01312022.rds')


###Markers from Anthony S
FeaturePlot(obj.seurat, features=c("PECAM1", "PDGFRA", "PDGFRB")) # Endothelial, Fibroblasts

#From Kai
FeaturePlot(obj.seurat, features=c("KRT4", "KRT15", "KRT14", "KRT10")) # epithelial markers

##Markers from Saba et al 2020 https://www-pnas-org.ezproxy.bu.edu/content/119/3/e2118424119/tab-figures-data
DefaultAssay(obj.seurat) <- 'RNA'
FeaturePlot(obj.seurat, features=c("CD74")) # Monocyte-derived macrophages and DC and LC
FeaturePlot(obj.seurat, features=c("MS4A7", "C1QA", "C1QB", "C1QC", "CCL5", "APOE")) # Monocyte-derived macrophages
FeaturePlot(obj.seurat, features=c("SIGLECH")) #pDC
FeaturePlot(obj.seurat, features=c("CCR7", "CCL22")) #migratory DCs
FeaturePlot(obj.seurat, features=c("EPCAM", "CD207")) # EPCAM+ DC, LC, respectively
FeaturePlot(obj.seurat, features=c("THY1", "CD3E")) # T cells
FeaturePlot(obj.seurat, features=c("NCR1")) # NK cells
FeaturePlot(obj.seurat, features=c("FOXP3", "CD8B1", "KLRK1", "TOP2A", "TRGC1")) # Tregs, CD8, activated cd4/cd8, proliferating, gamma delta; respectively


###########from source https://www.pnas.org/content/116/1/148
#Remarkably, the initial proposition of the involvement of EMT in metastasis was based on spatially varying levels of β-catenin (7), where cells at the invasive edge of the tumor and those in the tumor interior contained varying subcellular localization of β-catenin and E-cadherin, thus affecting cell adhesion and migration.
FeaturePlot(obj.seurat, features=c("CTNNB1", "CDH1")) # beta catenin gene, ecadherin
#the Notch-Jagged signaling pathway, which has been implicated in metastasis, drug resistance, and tumor relapse
#Notch-Jagged signaling among cancer cells can facilitate the formation of clusters of circulating tumor cells (CTCs)—the “bad actors” of cancer metastasis
#Notch-Jagged signaling between cancer cells and stromal cells and/or cells in metastatic niche can aggravate tumor progression
FeaturePlot(obj.seurat, features=c("NOTCH1", "JAG1"))
#many inflammatory cytokines such as IL-6, IL-1β, and TNF-α can promote Notch-Jagged signaling by increasing the intracellular production of Jagged
FeaturePlot(obj.seurat, features=c("IL6", "IL1B", "TNF"))
FeaturePlot(obj.seurat, features=c("SNAI1", "ZEB1")) #downstream of il6, up of notch and jag
#ALDH1+ epithelial-like BCSC population was localized in the interior close to the tumor stroma
FeaturePlot(obj.seurat, features=c("ALDH1A1"))

#to export
library(rBCS)
ExportSeurat(
  obj.seurat,
  "./sc4nqo02152022.bcs",
  unique.limit = 200, # Allow at most 200 unique in a categorical metadata
  clustering.name = "SCT_snn_res.0.8", # Use `my_cluster` as the clustering result
  batch.name = "treatment", # Use `sample_id` as batch identifiers
  compression.level = 9, # Maximize compression level
  author = "lkroeh@bu.edu", # I am the author
  raw.rna = "RNA", # Assay that has raw count of RNA
  norm.rna = "RNA", # Assay that has normalized count of RNA
  overwrite = TRUE # Overwite existing bcs file
)


#metadata
meta <- read.csv('./data/4nqo_annotation_and_metrics_12runs_combined.csv')

#subset control, 4nqo control, 4nqo+E7386 25mg, 4nqo+E7386 50mg, 4nqo+SP2509 (not drug controls)
#-2,3,4,10,11,12
Idents(obj.seurat) <- 'sample'
obj1 <- subset(x = obj.seurat, idents = c("KM_BM_0318_1C", "KM_BM_0318_5NC", "KM_BM_0318_6E", "KM_BM_0318_7E", "KM_BM_0318_8NL", "KM_BM_0319_1", "KM_BM_0319_5", "KM_BM_0319_6", "KM_BM_0319_7", "KM_BM_0319_8"), invert = FALSE)
Idents(obj1) <- 'seurat_clusters'
DimPlot(obj1, split.by = 'sample')
obj1$treatment <- factor(obj1$treatment, levels = c("control", "4NQO control","4NQO + E7386; 25 mg/kg", "4NQO+ E7386; 50 mg/kg", "4NQO+ SP2509"))
DimPlot(obj1, split.by = 'treatment')
#saveRDS(obj1, './data/allbutcontrol_02012022.rds')

#subset only control and cancer
Idents(obj.seurat) <- 'sample'
obj1 <- subset(x = obj.seurat, idents = c("KM_BM_0318_1C", "KM_BM_0318_5NC", "KM_BM_0319_1", "KM_BM_0319_5"), invert = FALSE)
obj1$treatment <- factor(obj1$treatment, levels = c("control", "4NQO control"))
#saveRDS(obj1, './data/Muzamilsubset_02092022.rds')

#subset only control and cancer and SP2509
Idents(obj.seurat) <- 'sample'
obj2 <- subset(x = obj.seurat, idents = c("KM_BM_0318_1C", "KM_BM_0318_5NC", "KM_BM_0319_1", "KM_BM_0319_5", "KM_BM_0319_8", "KM_BM_0318_8NL"), invert = FALSE)
obj2$treatment <- factor(obj2$treatment, levels = c("control", "4NQO control", "4NQO+ SP2509"))
#saveRDS(obj2, './data/02092022_SP2509subset.rds')
obj2 <- readRDS('./data/02092022_SP2509subset.rds')

obj2 <- RunPCA(obj2)
obj2 <- RunUMAP(obj2, dims = 1:30, verbose = FALSE)
obj2 <- FindNeighbors(obj2, dims = 1:30, verbose = FALSE)
obj2 <- FindClusters(obj2, verbose = FALSE)
DimPlot(obj2)

#subset only control and cancer and E7386
Idents(obj.seurat) <- 'sample'
obj2 <- subset(x = obj.seurat, idents = c("KM_BM_0318_1C", "KM_BM_0318_5NC", "KM_BM_0318_6E", "KM_BM_0318_7E", "KM_BM_0319_1", "KM_BM_0319_5", "KM_BM_0319_6", "KM_BM_0319_7"), invert = FALSE)
obj2$treatment <- factor(obj2$treatment, levels = c("control", "4NQO control","4NQO + E7386; 25 mg/kg", "4NQO+ E7386; 50 mg/kg"))

#prep for compass
Idents(obj1) <- 'singleRcelltype'
epi <- subset(x = obj1, idents = c("Epithelial cells"), invert = FALSE)
fib <- subset(x = obj1, idents = c("Fibroblasts"), invert = FALSE)
DimPlot(fib, split.by = 'treatment')

compassinput <- fib@assays$RNA@data
write.table(compassinput, file='./data/02072022_allsamples_FIBcellsonly_compassinput.tsv', quote=FALSE, sep='\t')

#Markers to send to Kai about epithelial cells
Idents(epi) <- 'seurat_clusters'
epimarkers <- FindAllMarkers(epi, min.pct = 0.5)
marks <- FindMarkers(epi, ident.1 = '10', ident.2 = '31')
write.csv(epimarkers, './results/epimarkers_02072022.csv')


##cell distribution plots by celltype
sce_df <- data.frame(obj2@meta.data[,c('sample', 'singleRcelltype', 'treatment')])
number_cells <- sce_df %>% 
  group_by(treatment) %>%
  summarise(number = n())
number_cells1 <- sce_df %>% 
  group_by(sample) %>%
  summarise(number = n())
DT::datatable(number_cells1)
grouped_labels2 <- sce_df %>% 
  group_by(singleRcelltype) %>%
  summarise(number = n())
DT::datatable(grouped_labels2)
grouped_labels1 <- sce_df %>% 
  group_by(singleRcelltype, treatment, sample) %>%
  summarise(number = n())
grouped_labels1$treatment <- factor(grouped_labels1$treatment, levels=c("control", "4NQO control","4NQO + E7386; 25 mg/kg", "4NQO+ E7386; 50 mg/kg", "4NQO+ SP2509"))


grouped_labels1 <- grouped_labels1 %>%
  mutate(numbernorm = case_when(
    (treatment == "4NQO control") ~ number/number_cells[which(number_cells$treatment == '4NQO control'),]$number,
    (treatment == "control") ~ number/number_cells[which(number_cells$treatment == 'control'),]$number,
    (treatment == "4NQO+ SP2509") ~ number/number_cells[which(number_cells$treatment == '4NQO+ SP2509'),]$number,
    (treatment == "4NQO + E7386; 25 mg/kg") ~ number/number_cells[which(number_cells$treatment == '4NQO + E7386; 25 mg/kg'),]$number,
    (treatment == "4NQO+ E7386; 50 mg/kg") ~ number/number_cells[which(number_cells$treatment == '4NQO+ E7386; 50 mg/kg'),]$number
  ))
grouped_labels1 <- grouped_labels1 %>%
  mutate(numbernorm = case_when(
    (sample == "KM_BM_0318_1C") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0318_1C'),]$number,
    (sample == "KM_BM_0318_5NC") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0318_5NC'),]$number,
    (sample == "KM_BM_0318_8NL") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0318_8NL'),]$number,
    (sample == "KM_BM_0319_1") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0319_1'),]$number,
    (sample == "KM_BM_0319_5") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0319_5'),]$number,
    (sample == "KM_BM_0319_8") ~ number/number_cells1[which(number_cells1$sample == 'KM_BM_0319_8'),]$number
  ))
grouped_labels3 <- grouped_labels1 %>%
  filter(singleRcelltype!='B cells, pro' & singleRcelltype!='Basophils' & singleRcelltype!='Microglia' & singleRcelltype!='NK cells')
ggplot2::ggplot(grouped_labels3, aes(singleRcelltype, numbernorm, color = treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="Cell-Types by Treatment(All)", y="# of Cells Normalized by Treatment",x="Cell-Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 8, face = "bold"))
ggplot2::ggplot(grouped_labels3, aes(singleRcelltype, log2(numbernorm), color = treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="Cell-Types by Treatment(All)", y="log2 # of Cells Normalized by Treatment",x="Cell-Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold")) #+ theme_bw()


##cell distribution plots by cluster
sce_df <- data.frame(obj1@meta.data[,c('sample', 'seurat_clusters', 'treatment')])
number_cells <- sce_df %>% 
  group_by(treatment) %>%
  summarise(number = n())
DT::datatable(number_cells)
grouped_labels2 <- sce_df %>% 
  group_by(seurat_clusters) %>%
  summarise(number = n())
DT::datatable(grouped_labels2)
grouped_labels1 <- sce_df %>% 
  group_by(seurat_clusters, treatment, sample) %>%
  summarise(number = n())
grouped_labels1$treatment <- factor(grouped_labels1$treatment, levels=c("control", "4NQO control","4NQO + E7386; 25 mg/kg", "4NQO+ E7386; 50 mg/kg", "4NQO+ SP2509"))

grouped_labels1 <- grouped_labels1 %>%
  mutate(numbernorm = case_when(
    (treatment == "4NQO control") ~ number/number_cells[which(number_cells$treatment == '4NQO control'),]$number,
    (treatment == "control") ~ number/number_cells[which(number_cells$treatment == 'control'),]$number,
    (treatment == "4NQO+ SP2509") ~ number/number_cells[which(number_cells$treatment == '4NQO+ SP2509'),]$number,
    (treatment == "4NQO + E7386; 25 mg/kg") ~ number/number_cells[which(number_cells$treatment == '4NQO + E7386; 25 mg/kg'),]$number,
    (treatment == "4NQO+ E7386; 50 mg/kg") ~ number/number_cells[which(number_cells$treatment == '4NQO+ E7386; 50 mg/kg'),]$number
  ))
ggplot2::ggplot(grouped_labels1, aes(seurat_clusters, numbernorm, color = treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="Cell-Types by Treatment(All)", y="# of Cells Normalized by Treatment",x="Cell-Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 8, face = "bold"))
ggplot2::ggplot(grouped_labels1, aes(seurat_clusters, log2(numbernorm), color = treatment)) +
  geom_boxplot(width=0.7) +
  labs(title="Cell-Types by Treatment(All)", y="# of Cells Normalized by Treatment",x="Cell-Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 8, face = "bold"))

##DE by treatment per cluster
Idents(obj1) <- 'seurat_clusters'
obj2 <- subset(x = obj1, idents = c("19"), invert = FALSE)
Idents(obj2) <- 'treatment'
DefaultAssay(obj2) <- 'RNA'
c12 <-FindAllMarkers(object =obj2, min.pct = 0.5, test.use = "MAST", slot = 'data') #, assay = sp@assays$originalexp@data, slot = 'data' for RNA

head(c12 %>%
       arrange(avg_log2FC, p_val_adj))
#DefaultAssay(g3) <- 'SCT'
FeaturePlot(obj1, features= c("FABP4"), split.by = 'treatment')
VlnPlot(obj1, features= c("SOCS3"), split.by = 'singleRcelltype', ncol = 2, pt.size = 0)

hm <- c12[(c12$avg_log2FC>=0.5 & c12$p_val_adj <0.05),]
DefaultAssay(obj2) <- 'SCT'
DoHeatmap(
  obj2,
  features = hm$gene,
  cells = NULL,
  group.by = "treatment",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = 'RNA',
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + theme(text = element_text(size = 7)) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

DT::datatable(hm)
