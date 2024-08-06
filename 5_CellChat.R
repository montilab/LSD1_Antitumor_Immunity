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
#saveRDS(cellchatem, './data/202408/cellchatcancer.rds')
#saveRDS(cellchatsp, './data/202408/cellchatcancerSP.rds')

cellchatem <- readRDS('./data/202408/cellchatcancer.rds')
cellchatsp <- readRDS('./data/202408/cellchatcancerSP.rds')

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
pdf('/restricted/projectnb/montilab-p/projects/oralcancer/4nqo/lina/results/202408/CellChat_CXCLcomparisonDC_toTCell.pdf', height = 4, width = 6)
netVisual_bubble(cellchat, sources.use = c("DCs", "ILC2"), targets.use = c("T/NK"), signaling = c("CXCL"), comparison = c(2, 1), angle.x = 45)
dev.off()
netVisual_bubble(cellchat, sources.use = c("DCs"), targets.use = c("T/NK"),  comparison = c(2, 1), max.dataset = 2, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)

#for circos
pos.dataset = "sp2509"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)

#net contains duplicates, we need to only include 1 when an interaction is common to both treatments
#example
net[net$source == "DCs" & net$target == "T/NK" & net$interaction_name == "CCL6_CCR2",]

#separate and combine 
netcancer <- net[net$datasets == 'cancer',]
netcancer$summarycol <- paste0(netcancer$source, netcancer$target, netcancer$interaction_name)
netsp <- net[net$datasets == 'sp2509',]
netsp$summarycol <- paste0(netsp$source, netsp$target, netsp$interaction_name)

#merge
netcomb <- merge(netcancer, netsp, by = "summarycol", all.x = T, all.y = T)
  
#add directionality
netcomb <- mutate(netcomb) %>%
  mutate(direction = case_when(
    (!is.na(prob.x) & !is.na(prob.y) & (prob.x > prob.y) ~ 2), #up in cancer
    (!is.na(prob.x) & !is.na(prob.y) & (prob.y > prob.x) ~ 1),  #up in SP
    (is.na(prob.x) & !is.na(prob.y) ~ 1), #only in SP
    (is.na(prob.y) & !is.na(prob.x) ~ 2) #only in cancer
)) 

#need consistent source name
netcomb <- mutate(netcomb) %>%
  mutate(sourcefinal = case_when(
    (is.na(source.x) ~ source.y),
    (is.na(source.y)  ~ source.x), 
    (!is.na(source.y) & !is.na(source.x)  ~ source.x)
)) 

#need consistent source name
netcomb <- mutate(netcomb) %>%
  mutate(targetfinal = case_when(
    (is.na(target.x) ~ target.y),
    (is.na(target.y)  ~ target.x),
    (!is.na(target.y) & !is.na(target.x)  ~ target.x)
  )) 

#get rid of cell types not in story( endo, fib, etc.)
netoi <- netcomb[(netcomb$sourcefinal == "Epithelial" | netcomb$sourcefinal == "Macrophages" | netcomb$sourcefinal == "DCs" | netcomb$sourcefinal == "ILC2" 
              | netcomb$sourcefinal == "Monocytes" | netcomb$sourcefinal == "T/NK" | netcomb$sourcefinal == "B Cells" | netcomb$sourcefinal == "Neutrophils") &
               (netcomb$targetfinal == "Epithelial" | netcomb$targetfinal == "Macrophages" | netcomb$targetfinal == "DCs" | netcomb$targetfinal == "ILC2" 
                | netcomb$targetfinal == "Monocytes" | netcomb$targetfinal == "T/NK" | netcomb$targetfinal == "B Cells" | netcomb$targetfinal == "Neutrophils"), ]

FC <- netoi
FC$sourcefinal <- factor(FC$sourcefinal, levels = c(unique(FC$sourcefinal)))
FC$targetfinal <- factor(FC$targetfinal, levels = c(unique(FC$targetfinal)))

#FC$direction <- factor(FC$direction, levels = c(1,2))

links <- data.frame(from = FC$sourcefinal,
                    to = FC$targetfinal,
                    value = as.numeric(FC$direction))

library(circlize)
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

write.csv(netcomb, './results/202408/cellchat_comparison.csv')

