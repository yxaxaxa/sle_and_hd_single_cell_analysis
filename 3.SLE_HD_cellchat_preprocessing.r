library(Seurat)
library(ggplot2)
library(Matrix)

#SLE flare data read in
#can be accessed through https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137029
mm<-readMM('./project/sle_flare_GSE137029/matrix.mtx.gz')
features=read.table('./project/sle_flare_GSE137029/features.tsv.gz',sep='\t')
bar<-read.table('./project/sle_flare_GSE137029/barcodes.tsv.gz',sep='\t')
rownames(mm)<-bar$V1
colnames(mm)<-features$V1
sle<-CreateSeuratObject(counts=t(mm),project='sle')

#HD Sample 
##can be accessed through https://cells.ucsc.edu/?ds=multimodal-pbmc+sct data download page
##hd rds
##first time processing hd data please run following step for first hd data processing
##mat <- fread("exprMatrix.tsv.gz")
##meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
##genes = mat[,1][[1]]
##genes = gsub(".+[|]", "", genes)
##mat = data.frame(mat[,-1], row.names=genes)
##hd <- CreateSeuratObject(counts = mat, project = "hd", meta.data=meta)
##hd@meta.data$orig.ident<-'hd'
##saveRDS(hd,ucsc_pbmc_hd.rds)
hd<-readRDS('./project/ucsc_pbmc_hd.rds')


####child sle used datasets
path<-'./project/GSE135779_SLE/cSLE/SLE_HEAVY/'
file<-dir(path)
file

file.path<-lapply(file,function(x){
    paste(path,x,sep='')
})
csle_heavy<-list()
for(i in 1:length(file.path)){
    csle_heavy[[i]]<-Read10X(file.path[[i]])
    csle_heavy[[i]]<-CreateSeuratObject(csle_heavy[[i]],project='csle_heavy')
    csle_heavy[[i]]@meta.data$orig.ident<-'csle_heavy'
    
    
}


csle_heavy_m<-csle_heavy[[1]]
for(i in 2:length(csle_heavy)){
    csle_heavy_m<-merge(csle_heavy_m,csle_heavy[[i]])
}

sle<-merge(sle,csle_heavy_m)

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
data.integrated<-FindClusters(data.integrated,res=0.2)

DimPlot(data.integrated,group.by=c('seurat_clusters','celltype.l1'),label=T,raster=F)

new.cluster.ids <- c("T", "Mono", "T", "T", "B", "NK",
    "Mono", "Mono", "Mono",'DC','other','DC','other','T','other','T','other','DC','B','other')
names(new.cluster.ids) <- levels(data.integrated)
data.integrated <- RenameIdents(data.integrated, new.cluster.ids)
DimPlot(data.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(data.integrated,'./project/sle_nature_communiation_and_immunnology_merge_include_hd.rds')
