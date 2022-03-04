library(Seurat)
library(ggplot2)
library(Matrix)

#SLE data import
mm<-readMM('./project/sle_flare_GSE137029/matrix.mtx.gz')
features=read.table('./project/sle_flare_GSE137029/features.tsv.gz',sep='\t')
bar<-read.table('./project/sle_flare_GSE137029/barcodes.tsv.gz',sep='\t')
rownames(mm)<-bar$V1
colnames(mm)<-features$V1
sle<-CreateSeuratObject(counts=t(mm),project='sle')

###HD Sample 
###can be accessed through https://cells.ucsc.edu/?ds=multimodal-pbmc+sct data download page
###for exprsmatrix:wget https://cells.ucsc.edu/multimodal-pbmc/sct/exprMatrix.tsv.gz
###for metadata :wget https://cells.ucsc.edu/multimodal-pbmc/sct/meta.tsv
###hd rds
###first time processing hd data please run following step for first hd data processing
###mat <- fread("exprMatrix.tsv.gz")
###meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
###genes = mat[,1][[1]]
###genes = gsub(".+[|]", "", genes)
###mat = data.frame(mat[,-1], row.names=genes)
###hd <- CreateSeuratObject(counts = mat, project = "hd", meta.data=meta)
###hd@meta.data$orig.ident<-'hd'
###saveRDS(hd,ucsc_pbmc_hd.rds)
hd<-readRDS('./project/ucsc_pbmc_hd.rds')

### data merge with RPCA
data.list<-list(sle,hd)
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = data.list, reference = c(1,2), reduction = "rpca",
    dims = 1:50)
data.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, dims = 1:50)
data.integrated<-FindNeighbors(data.integrated,dims=1:50)


DefaultAssay(data.integrated)<-'integrated'
data.integrated<-FindClusters(data.integrated,res=0.1)

orig<-data.integrated@meta.data$orig.ident
orig[orig!='sle']<-'HD'
orig[orig=='sle']<-'SLE'
data.integrated@meta.data$orig.ident<-orig

p1<-DimPlot(data.integrated,group.by=c('orig.ident'),raster=F,label=T)
p1

p2<-DimPlot(data.integrated,group.by=c("seurat_clusters"))
p2

### save merged data as RDS file
saveRDS(data.integrated,'./project/sle_hd_merge_1126.rds')
