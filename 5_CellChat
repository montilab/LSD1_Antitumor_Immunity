library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(patchwork)

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina')
obj.seurat <- readRDS('./data/202408/20240801_SP2509.rds')
imm <- readRDS('./data/202408/imm.Rds')

#set up labels
obj.seurat$ccclabels1 <- obj.seurat$CellType
obj.seurat$ccclabels2 <- imm$ImmIdents

df <- obj.seurat@meta.data
df <- df %>% mutate(ccclabelsx =
                      case_when(ccclabels2 == "TGD" ~ "T/NK",
                                ccclabels2 != "TGD" ~ ccclabels2,
                                ccclabels1 == "Epithelial" ~ ccclabels1,
                                ccclabels1 == "Endothelial" ~ ccclabels1,
                                ccclabels1 == "Lymphatic Endothelial" ~ ccclabels1,
                                ccclabels1 == "Fibroblasts" ~ ccclabels1,
                                ccclabels1 == "StromalVSMC" ~ ccclabels1,
                                ccclabels1 == "Glial" ~ ccclabels1,
                                ccclabels1 == "RBCs" ~ ccclabels1))
                                #is.na(ccclabels2) ~ ccclabels1))

obj.seurat <- AddMetaData(obj.seurat, df$ccclabelsx, col.name = 'ccclabels')
DimPlot(obj.seurat, group.by = 'ccclabels', label = T)
Idents(obj.seurat) <- 'ccclabels'

#split treatments
Idents(obj.seurat) <- 'treatment'
nqoctrl <- subset(obj.seurat, idents = c("4NQO control"))
sptrt <- subset(obj.seurat, idents = c("4NQO+ SP2509"))

#add sample info
nqoctrl$samples <- nqoctrl$sample
sptrt$samples <- sptrt$sample

library(stringr)
#need to make gene names not uppercase
nqoctrl@assays$RNA@counts@Dimnames[[1]] <- str_to_title(nqoctrl@assays$RNA@counts@Dimnames[[1]])
nqoctrl@assays$RNA@data@Dimnames[[1]] <- str_to_title(nqoctrl@assays$RNA@data@Dimnames[[1]])

sptrt@assays$RNA@counts@Dimnames[[1]] <- str_to_title(sptrt@assays$RNA@counts@Dimnames[[1]])
sptrt@assays$RNA@data@Dimnames[[1]] <- str_to_title(sptrt@assays$RNA@data@Dimnames[[1]])

###run
#control cellchat
CellChatDB <- CellChatDB.mouse
cellchatem <- createCellChat(nqoctrl, group.by = 'ccclabels', assay = 'RNA')
cellchatem <- setIdent(cellchatem, ident.use = "ccclabels")
levels(cellchatem@idents)
groupSizeem <- as.numeric(table(cellchatem@idents))

#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB
cellchatem@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchatem <- subsetData(cellchatem) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 1) # do parallel
cellchatem <- identifyOverExpressedGenes(cellchatem, thresh.pc = 0.01, thresh.fc = 0.01)
cellchatem <- identifyOverExpressedInteractions(cellchatem)

#cellchatem <- computeCommunProb(cellchatem,type = "truncatedMean", trim = 0)
cellchatem <- computeCommunProb(cellchatem, trim = 0.05, type = "truncatedMean")
#computeAveExpr(cellchatem, features = c("Csf1","Il34"), type =  "truncatedMean", trim = 0.1)
cellchatem <- filterCommunication(cellchatem, min.cells = 5)

df.netem <- subsetCommunication(cellchatem)

cellchatem <- computeCommunProbPathway(cellchatem)

cellchatem <- aggregateNet(cellchatem)

groupSizeem <- as.numeric(table(cellchatem@idents))

cellchatem <- netAnalysis_computeCentrality(cellchatem)

##sp cellchat
CellChatDB <- CellChatDB.mouse
cellchatsp <- createCellChat(sptrt, group.by = 'ccclabels', assay = 'RNA')
cellchatsp <- setIdent(cellchatsp, ident.use = "ccclabels")
levels(cellchatsp@idents)
groupSizesp <- as.numeric(table(cellchatsp@idents))

#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB
cellchatsp@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchatsp <- subsetData(cellchatsp) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel
cellchatsp <- identifyOverExpressedGenes(cellchatsp, thresh.pc = 0.01, thresh.fc = 0.01)
cellchatsp <- identifyOverExpressedInteractions(cellchatsp)


#cellchatsp <- computeCommunProb(cellchatsp,type = "truncatedMean", trim = 0)
cellchatsp <- computeCommunProb(cellchatsp, trim = 0.05, type = "truncatedMean")#, population.size = FALSE, raw.use = FALSE, nboot = 1)
#cellchatsp@net
#computeAveExpr(cellchatsp, features = c("Csf1","Il34"), type =  "truncatedMean", trim = 0.1)
cellchatsp <- filterCommunication(cellchatsp, min.cells = 5)

df.netsp <- subsetCommunication(cellchatsp)

cellchatsp <- computeCommunProbPathway(cellchatsp)

cellchatsp <- aggregateNet(cellchatsp)

groupSizesp <- as.numeric(table(cellchatsp@idents))

cellchatsp <- netAnalysis_computeCentrality(cellchatsp)


###done 
saveRDS(cellchatem, './data/202408/cellchatcancer.rds')
saveRDS(cellchatsp, './data/202408/cellchatcancerSP.rds')

#check for specific interactions
cellchatem@net$prob[, , "CXCL9_CXCR3"]["DCs", "T/NK"]
cellchatsp@net$prob[, , "CXCL9_CXCR3"]["DCs", "T/NK"]

cellchatem@net$pval[, , "CXCL9_CXCR3"]["DCs", "T/NK"]
cellchatsp@net$pval[, , "CXCL9_CXCR3"]["DCs", "T/NK"]


#combine
object.list <- list(cancer = cellchatem, sp2509 = cellchatsp)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#individual
netVisual_bubble(cellchatem, sources.use = c("DCs"), targets.use = c("T/NK"), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchatsp, sources.use = c("DCs"), targets.use = c("T/NK"), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#comparison
netVisual_bubble(cellchat, sources.use = c("DCs", "ILC2"), targets.use = c("T/NK"), signaling = c("CXCL"), comparison = c(2, 1), angle.x = 45)

#for circos
pos.dataset = "sp2509"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)

#get rid of cell types not in story( endo, fib, etc.)
netoi <- net[(net$source == "Epithelial" | net$source == "Macrophages" | net$source == "DCs" | net$source == "ILC2" 
              | net$source == "Monocytes" | net$source == "T/NK" | net$source == "B Cells" | net$source == "Neutrophils") &
               (net$target == "Epithelial" | net$target == "Macrophages" | net$target == "DCs" | net$target == "ILC2" 
                | net$target == "Monocytes" | net$target == "T/NK" | net$target == "B Cells" | net$target == "Neutrophils"), ]

FC <- netoi
FC$source <- factor(FC$source, levels = c(unique(FC$source)))
FC$target <- factor(FC$target, levels = c(unique(FC$target)))

#add direction based on up in cancer or up in sp2509
FC <- mutate(FC) %>%
  mutate(direction = case_when(
    (datasets == "sp2509") ~ 1, #up in sp
    (datasets == "cancer") ~ 2)) #up in cancer

FC$direction <- factor(FC$direction, levels = c(1,2))

links <- data.frame(from = FC$source,
                    to = FC$target,
                    value = as.numeric(FC$direction))

col_fun = colorRamp2(c(1,2), c("blue", "red"), transparency = 0.5)

chordDiagram(links, directional = T, self.link = 1,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", col = col_fun)#,
             #grid.col = c( "hotpink","seagreen3", "steelblue1","yellowgreen","pink2", "gold2", "orange", "red2") ) #link.target.prop = FALSE , scale = T



pdf('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina/results/202408/circosAllCells.pdf', height = 4, width = 4)
chordDiagram(links, directional = T, self.link = 1,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", col = col_fun)#,
#grid.col = c( "hotpink","seagreen3", "steelblue1","yellowgreen","pink2", "gold2", "orange", "red2") ) #link.target.prop = FALSE , scale = T
dev.off()

pdf('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina/results/202408/circosImmCells.pdf', height = 4, width = 4)
chordDiagram(links, directional = T, self.link = 1,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", col = col_fun,
             grid.col = c( "hotpink","seagreen3", "steelblue1","yellowgreen","pink2", "gold2", "orange", "red2") ) #link.target.prop = FALSE , scale = T
dev.off()

write.csv(net, './results/202408/cellchat_comparison.csv')

