library(Seurat)
library(CellChat)
library(pheatmap)

sle<-readRDS('./project/sle_hd_merge_1126.rds')

options(repr.plot.width=8,repr.plot.height=8)
sle<-FindClusters(sle,res=0.3)
DimPlot(sle,group.by=c('seurat_clusters'),label=T,raster=F,label.size=6,repel=T,pt.size=0.01)

###b
DefaultAssay(sle)<-'integrated'
sle_sub_B<-subset(sle,seurat_clusters%in%c(6,10))
sle_sub_B<-RunPCA(sle_sub_B)
sle_sub_B <- RunUMAP(object = sle_sub_B, dims = 1:20)
sle_sub_B <- FindNeighbors(sle_sub_B , reduction = "pca", dims = 1:20)
sle_sub_B<-FindClusters(sle_sub_B,res=0.3)

#dir.create("./sle_picture/")
#dir.create("./sle_picture/B/")
options(repr.plot.width=17.6,repr.plot.height=8)
sle_sub_B@meta.data$orig.ident<-toupper(sle_sub_B@meta.data$orig.ident)
p1<-DimPlot(sle_sub_B,group.by='seurat_clusters',label=F,raster=F,label.size=6,repel=T,pt.size=0.01,cols=c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9',
                                            '#2f7fc1','#82b1d2','#96c37d','#f4d267'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle_sub_B,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.01,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p<-p1+p2
p
#options(repr.plot.width=17.6,repr.plot.height=8)
#p1
#p2
#ggsave(p, file="./sle_picture/B/B_dimplot.pdf",width = 17.6, height = 8)
#ggsave(p1, file="./sle_picture/B/B_p1_dimplot.tif",device = "tiff",width = 5, height =5)
#ggsave(p2, file="./sle_picture/B/B_p2_dimplot.tif",device = "tiff",width = 5, height =5)


options(repr.plot.width=17.6,repr.plot.height=8)
count<-as.data.frame(table(sle_sub_B@meta.data$seurat_clusters,sle_sub_B@meta.data$orig.ident))
count$Var2<-toupper(count$Var2)
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count)+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,size=3,color='white',lty=2)+
theme( 
      axis.text.x=element_text(size=45),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=45),
     axis.text.y=element_text(size=45),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))+labs(title='')&NoLegend()
p
#ggsave(p, file="./sle_picture/B/B_bar.pdf",width = 17.6, height = 6)

options(repr.plot.width=25,repr.plot.height=15)
DefaultAssay(sle_sub_B)<-'RNA'
features=c('CD19',"MS4A1","ISG15","IFI44L","IFI27","MX1","IGHD","CD27","CD24","CD38","TBX21","ITGAX","CXCR5","TRAF5","CR2")
p<-FeaturePlot(sle_sub_B,features=features,cols=c('#eee6e3','#bd3c39'),pt.size=.5,ncol=5)&NoAxes()&NoLegend()
p
#ggsave(p, file="./sle_picture/B/b_feature.pdf",width =25, height = 15)
#ggsave(p, file="./sle_picture/B/b_feature.tif",device = "tiff",width = 25, height = 15)

markerb<-FindAllMarkers(sle_sub_B,only.pos=T)
#write.csv(markerb,'./sle_b_cell_cluster_marker.csv')

options(repr.plot.width=7.93,repr.plot.height=9.47)
DefaultAssay(sle_sub_B)<-'integrated'
markerb<-read.csv('./sle_b_cell_cluster_marker.csv')
genea<-markerb%>%  group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
p<-DoHeatmap(subset(sle_sub_B,downsample=300),features=genea$gene,group.colors =c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9', '#2f7fc1','#82b1d2','#96c37d','#f4d267'),angle = 0)+scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")) ) + guides(color="none")
p
#ggsave(p, file="./sle_picture/B/B_heatmap.pdf",width = 7.93, height = 9.47)

###mono
DefaultAssay(sle)<-'integrated'
sle_sub_Mono<-subset(sle,seurat_clusters%in%c(0,7,8,9))
sle_sub_Mono<-RunPCA(sle_sub_Mono)
sle_sub_Mono <- RunUMAP(object = sle_sub_Mono, dims = 1:20)
sle_sub_Mono <- FindNeighbors(sle_sub_Mono , reduction = "pca", dims = 1:20)
sle_sub_Mono<-FindClusters(sle_sub_Mono,res=0.3)

#dir.create("./sle_picture/Mono/")
options(repr.plot.width=17.6,repr.plot.height=8)
DefaultAssay(sle_sub_Mono)<-'integrated'
sle_sub_Mono<-subset(sle_sub_Mono,seurat_clusters!=9)
sle_sub_Mono@meta.data$orig.ident<-toupper(sle_sub_Mono@meta.data$orig.ident)
p1<-DimPlot(sle_sub_Mono,group.by='seurat_clusters',label=F,raster=F,label.size=6,repel=T,pt.size=0.01,cols=c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9',
                                            '#2f7fc1','#82b1d2','#96c37d','#f4d267'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle_sub_Mono,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.01,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')

p <- p1+p2
p
#&NoLegend()
#options(repr.plot.width=5,repr.plot.height=5)
#p1
#p2
#ggsave(p, file="./sle_picture/B/B_dimplot.pdf",width = 17.6, height = 8)
#ggsave(p1, file="./sle_picture/Mono/Mono_p1_dimplot.tif",device = "tiff",width = 5, height =5)
#ggsave(p2, file="./sle_picture/Mono/Mono_p2_dimplot.tif",device = "tiff",width = 5, height =5)

options(repr.plot.width=17.6,repr.plot.height=8)
count<-as.data.frame(table(sle_sub_Mono@meta.data$seurat_clusters,sle_sub_Mono@meta.data$orig.ident))
count$Var2<-toupper(count$Var2)
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count[count$Var1!='9',])+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,size=3,color='white',lty=2)+
theme( 
      axis.text.x=element_text(size=45),
      axis.title.x=element_blank(),
     # axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=45),
     axis.text.y=element_text(size=45),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))+labs(title='')&NoLegend()
p
#ggsave(p,file="./sle_picture/Mono/Mono_bar.pdf",width = 17.6, height = 6)

marker_Mono<-FindAllMarkers(sle_sub_Mono,only.pos=T)
#write.csv(marker_Mono,'./differential_marker_Mono.csv')

options(repr.plot.width=7.93,repr.plot.height=8.74)




DefaultAssay(sle_sub_Mono)<-'integrated'
marker_Mono<-read.csv('./differential_marker_mono.csv',row.names=1)
geneb<-marker_Mono%>%  group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
p<-DoHeatmap(subset(sle_sub_Mono,downsample=300),features=geneb$gene,group.colors =c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9','#2f7fc1','#82b1d2','#96c37d','#f4d267'),angle=0)+scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")) ) + guides(color="none")
p
#ggsave(p,file="./sle_picture/Mono/Mono_heatmap.pdf",width = 8.67, height = 9.3)

options(repr.plot.width=25,repr.plot.height=15)
DefaultAssay(sle_sub_Mono)<-'RNA'
features=c('LYZ','S100A8','CD14','FCGR3A','ISG15','IFI27','IFI44L','MX2','FKBP5')
p<-FeaturePlot(sle_sub_Mono,features=features,cols=c('#eee6e3','#bd3c39'),pt.size=.5,ncol=4)&NoLegend()&NoAxes()
p
#ggsave(p,file="./sle_picture/Mono/Mono_feature.pdf",width = 25, height = 15)
#ggsave(p,file="./sle_picture/Mono/Mono_feature.tif",device="tiff",width = 25, height = 15)

####t
DefaultAssay(sle)<-'integrated'
sle_sub<-subset(sle,seurat_clusters%in%c(1,2,3,5,11,14))
sle_sub<-RunPCA(sle_sub)
sle_sub <- RunUMAP(object = sle_sub, dims = 1:20)
sle_sub <- FindNeighbors(sle_sub , reduction = "pca", dims = 1:20)
sle_sub<-FindClusters(sle_sub,res=0.25)

options(repr.plot.width=17.6,repr.plot.height=8)
sle_sub@meta.data$orig.ident<-toupper(sle_sub@meta.data$orig.ident)
#dir.create("./sle_picture/T/")
p1<-DimPlot(sle_sub,group.by='seurat_clusters',label=F,raster=F,label.size=6,repel=T,pt.size=0.01,cols=c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9',
                                            '#2f7fc1','#82b1d2','#96c37d','#f4d267','#e88482','#c498b2', '#9394e7' ))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle_sub,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.01,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p<-p1+p2
p
#options(repr.plot.width=5,repr.plot.height=5)
#p1
#p2
#ggsave(p1, file="./sle_picture/T/T_p1_dimplot.tiff",device="tiff",width = 5, height = 5)
#ggsave(p2, file="./sle_picture/T/T_p2_dimplot.tiff",device="tiff",width = 5, height = 5)

options(repr.plot.width=8,repr.plot.height=5)
count<-as.data.frame(table(sle_sub@meta.data$seurat_clusters,sle_sub@meta.data$orig.ident))
count$Var2<-toupper(count$Var2)
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count)+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,color='white',size=3,lty=2)+
theme( 
      axis.text.x=element_text(size=20),
      axis.title.x=element_blank(),
     # axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=20),
     axis.text.y=element_text(size=20),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))
p
#ggsave(p, file="./sle_picture/T/T_bar.pdf",width = 8, height = 5)

####nk
DefaultAssay(sle)<-'integrated'
sle_sub1<-subset(sle,seurat_clusters%in%c(4,17))
sle_sub1<-RunPCA(sle_sub1)
sle_sub1 <- RunUMAP(object = sle_sub1, dims = 1:20)
sle_sub1 <- FindNeighbors(sle_sub1 , reduction = "pca", dims = 1:20)
sle_sub1<-FindClusters(sle_sub1,res=0.1)

options(repr.plot.width=17.6,repr.plot.height=8)
dir.create("./sle_picture/NK/")
sle_sub1@meta.data$orig.ident<-toupper(sle_sub1@meta.data$orig.ident)
p1<-DimPlot(sle_sub1,group.by='seurat_clusters',label=F,raster=F,label.size=6,repel=T,pt.size=0.01,cols=c('#f7e1ed','#b28247','#d86363','#44a064','#8dcfc9',
                                            '#2f7fc1','#82b1d2','#96c37d','#f4d267'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle_sub1,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.01,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')

p<-p1+p2
p
#options(repr.plot.width=5,repr.plot.height=5)
#p1
#p2
ggsave(p1, file="./sle_picture/NK/NK_dimplot.tif",device="tiff",width = 5, height = 5)
ggsave(p2, file="./sle_picture/NK/NK_dimplot.tif",device="tiff",width = 5, height = 5)

options(repr.plot.width=8,repr.plot.height=5)
count<-as.data.frame(table(sle_sub1@meta.data$seurat_clusters,sle_sub1@meta.data$orig.ident))
count$Var2<-toupper(count$Var2)
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count)+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,color='white',size=3,lty=2)+
theme( 
      axis.text.x=element_text(size=20),
      axis.title.x=element_blank(),
     # axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=20),
     axis.text.y=element_text(size=20),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))
p
#ggsave(p, file="./sle_picture/NK/NK_bar.pdf",width = 8, height = 5)

####dc
DefaultAssay(sle)<-'integrated'
sle_sub4<-subset(sle,seurat_clusters%in%c(12,16))
sle_sub4<-RunPCA(sle_sub4)
sle_sub4 <- RunUMAP(object = sle_sub4, dims = 1:20)
sle_sub4 <- FindNeighbors(sle_sub4 , reduction = "pca", dims = 1:20)
sle_sub4<-FindClusters(sle_sub4,res=0.1)

options(repr.plot.width=17.6,repr.plot.height=8)
#dir.create("./sle_picture/DC/")
sle_sub4@meta.data$orig.ident<-toupper(sle_sub4@meta.data$orig.ident)
p1<-DimPlot(sle_sub4,group.by='seurat_clusters',label=F,raster=F,label.size=6,repel=T,pt.size=0.1,cols=c('#f7e1ed','#b28247','#d86363'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle_sub4,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.1,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p <- p1+p2
p
#options(repr.plot.width=5,repr.plot.height=5)
#p1
#p2
#ggsave(p1,file="./sle_picture/DC/DC_p1_heatmap.tif",width = 5, height = 5)
#ggsave(p2,file="./sle_picture/DC/DC_p2_heatmap.tif",width = 5, height = 5)

options(repr.plot.width=4,repr.plot.height=5)
count<-as.data.frame(table(sle_sub4@meta.data$seurat_clusters,sle_sub4@meta.data$orig.ident))
count$Var2<-toupper(count$Var2)
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count)+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,color='white',size=3,lty=2)+
theme( 
      axis.text.x=element_text(size=20),
      axis.title.x=element_blank(),
     # axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=20),
     axis.text.y=element_text(size=20),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))
p
#ggsave(p,file="./sle_picture/DC/DC_bar.pdf",width = 4, height = 5)

options(repr.plot.width=17.6,repr.plot.height=8)
#dir.create("./sle_picture/ALL/")
sle<-FindClusters(sle,res=0.1)
orig<-sle@meta.data$orig.ident
orig[orig=='sle']<-'SLE'
orig[orig=='hd']<-'HD'
sle@meta.data$orig.ident<-orig
p1<-DimPlot(sle, ,label=F,raster=F,label.size=6,repel=T,pt.size=0.01,cols=c("#f7e1ed",'#b28247','#d86363','#44a064','#8dcfc9',
                                            '#2f7fc1','#82b1d2','#96c37d','#f4d267','#e88482','#c498b2', '#9394e7' ,'#e7dad2','#fea3a2','#dae06d','#5c85a5'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')
p2<-DimPlot(sle,group.by=c('orig.ident'),label=F,raster=F,label.size=4,repel=T,pt.size=0.01,cols=c('#d86363','#2f7fc1'))+theme(
                                                                              axis.title.x=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.ticks.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks.y=element_blank(),
                                                                              axis.line=element_blank())+labs(title='')



p <- p1+p2
p
#options(repr.plot.width=5,repr.plot.height=5)
#p1
#p2
#ggsave(p1,file="./sle_picture/ALL/all_p1_dimplot.tif",device="tiff",width = 5, height = 5)
#ggsave(p2,file="./sle_picture/ALL/all_p2_dimplot.tif",device="tiff",width = 5, height = 5)


options(repr.plot.width=16,repr.plot.height=5)
count<-as.data.frame(table(sle@meta.data$seurat_clusters,sle@meta.data$orig.ident))
colnames(count)<-c('Var1','Category','Freq')
p<-ggplot(count)+geom_bar(aes(x=Var1,y=Freq,fill=Category),stat='identity',position = 'fill')+theme_test()+geom_hline(yintercept=0.75,size=3,color='white',lty=2)+
theme( 
      axis.text.x=element_text(size=20),
      axis.title.x=element_blank(),
     # axis.ticks.x=element_blank(),
      axis.title.y=element_text(size=20),
     axis.text.y=element_text(size=15),
    legend.title=element_text(size=15),
legend.text=element_text(size=12),legend.position = 'right')+
scale_fill_manual(values=c('#d86363','#2f7fc1'))
p
#ggsave(p,file="./sle_picture/ALL/all_bar.pdf",width = 16, height = 5)

options(repr.plot.width=17.6,repr.plot.height=7)
gene<-c('HBB','PPBP','S100A8','CD14','LYZ','CST3','CD1C','FCER1G','FCGR3A','IL7R','CD3E','CD3D','CD8A','GZMA','GZMB',
       'GZMH','NKG7','MS4A1','CD79A','MZB1','LILRA4')
p<-DotPlot(sle,features=gene,cols=c('darkgrey','red'),scale = T)
p
ggsave(p,file="./sle_picture/ALL/all_dotplot.pdf",width = 16, height = 7)
