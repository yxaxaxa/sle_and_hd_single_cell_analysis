library(dplyr)
marker_mono<-read.csv('./differential_marker_mono',row.names=1)
geneb<-marker_mono%>%  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
markerb<-read.csv('./sle_b_cell_cluster_marker.csv')
genea<-markerb%>%  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

length(c(unique(genea$gene),unique(geneb$gene),'ICAM1','ICAM2','ITGAL','SPN','ITGAX','ITGAM','APP','CD74','THBS1','CD47','CD36',
         'IL16','TNFRSF1A','TNF','TGFB1','TGFBR2','SORT1','BAG6'))
####combine B cell cluster marker monocyte cluster marker and Cellchat identified marker to a total 86 gene set
write.csv(unique(c(unique(genea$gene),unique(geneb$gene),'ICAM1','ICAM2','ITGAL','SPN','ITGAX','ITGAM','APP','CD74','THBS1','CD47','CD36',
                   'IL16','TNFRSF1A','TNF','TGFB1','TGFBR2','SORT1','BAG6')),'three_gene_combined_markerset.csv')