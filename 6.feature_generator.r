library(dplyr)
library(stringr)
library(readxl)

markerb<-read.csv('../script4paper/differential_marker_b.csv')
marker_mono<-read.csv('../script4paper/differential_marker_mono.csv')
genea<-markerb%>%  group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
geneb<-marker_mono%>%  group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)



geneset2<-read_excel('./CELLCHAT_GENE_hd_specific.xlsx')

gene_int<-unique(c(genea$gene,geneb$gene))

genec<-geneset2[geneset2$count>=3,]$`...1`

length(c(gene_int,genec))
length(unique(gene_int,genec))

gene<-c(genea$gene,geneb$gene,genec)

write.csv(unique(genec),'./combined_gene_for_machine_learning_cellchat.csv')


