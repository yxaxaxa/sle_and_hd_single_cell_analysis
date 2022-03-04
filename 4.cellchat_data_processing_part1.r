library(Seurat)
library(CellChat)
library(pheatmap)


sle<-readRDS('./project/sle_nature_communiation_and_immunnology_merge_include_hd.rds')

orig<-sle@meta.data$orig.ident
orig[orig!='csle_heavy'&orig!='sle']<-'hd'
sle@meta.data$orig.ident<-orig
table(sle@meta.data$orig.ident)

hd_group<-subset(sle,orig.ident=='hd')

labels <- Idents(hd_group)
hd_data.input <- GetAssayData(hd_group, assay = "RNA", slot = "data") # normalized data matrix
#labels <- celltype
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
hd_cellchat <- createCellChat(object = hd_data.input, meta = meta, group.by = "group")
hd_cellchatDB <- CellChatDB.human  # use hd_cellchatDB.mouse if running on mouse data
showDatabaseCategory(hd_cellchatDB)
hd_cellchatDB.use <- hd_cellchatDB
hd_cellchat@DB <- hd_cellchatDB.use
hd_cellchat <- subsetData(hd_cellchat) 
future::plan("multicore", workers = 4)
hd_cellchat <- identifyOverExpressedGenes(hd_cellchat)
hd_cellchat <- identifyOverExpressedInteractions(hd_cellchat)
hd_cellchat <- projectData(hd_cellchat, PPI.human)
options(future.globals.maxSize= 1891289600) 
hd_cellchat <- computeCommunProb(hd_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
hd_cellchat <- filterCommunication(hd_cellchat, min.cells = 10)
hd_cellchat <- computeCommunProbPathway(hd_cellchat)
hd_cellchat <- aggregateNet(hd_cellchat)
hd_cellchat <- netAnalysis_computeCentrality(hd_cellchat, slot.name = "netP")

sle_group<-subset(sle,orig.ident=='sle')

#sle_group<-sle_group[,sample(1:dim(sle_group)[2],150000)]
data.input <- GetAssayData(sle_group, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(sle_group)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
options(future.globals.maxSize= 1891289600) 
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
sle_cellchat<-cellchat

csle_heavy_group<-subset(sle,orig.ident=='csle_heavy')

#csle_heavy_group<-csle_heavy_group[,sample(1:dim(csle_heavy_group)[2],150000)]
data.input <- GetAssayData(csle_heavy_group, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(csle_heavy_group)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multicore", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
options(future.globals.maxSize= 1891289600) 
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
csle_heavy_cellchat<-cellchat

object.list<-list(hd_cellchat,sle_cellchat,csle_heavy_cellchat)
saveRDS(object.list,'./project/merge_sle_immunonolgy_communication_cellchat.rds')


