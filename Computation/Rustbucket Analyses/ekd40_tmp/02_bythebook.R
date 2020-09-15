#!/usr/bin/env R

library(dplyr)
library(tibble)
library(Seurat); s_o <- readRDS('s_o.clustered.RDS')
s_o <- RunTSNE(s_o, dims.use = 1:45, max_iter=20000, perplexity=5, check_dupli
		 cates = FALSE, theta=0.05)
s_o <- FindNeighbors(s_o, dims = 1:45)
s_o <- FindClusters(s_o, resolution = 12)
pdf('tsne.pdf', width=8, height=4)
DimPlot(s_o, reduction = "tsne", label=T)
dev.off()
pdf('tsne-features.pdf', width=8, height=4)
FeaturePlot(s_o, features = c("SoxN", "D", "vnd", "wor", "pros"))
dev.off()
saveRDS(s_o, 's_o.bythebook.RDS')


s_o.markers <- FindAllMarkers(s_o, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(s_o.markers, 'so_markers.bythebook.tsv', sep='\t')
write.csv(s_o.markers, 'ClusterMap/allen.markers.csv')

pdf('violin_bythebook.pdf')
RidgePlot(s_o, features = c("SoxN", "D"), ncol = 2)
dev.off()



s_o <- RunUMAP(s_o, dims = 1:45)
pdf('umap.pdf', width=8, height=4)
DimPlot(s_o, reduction = 'umap', label=T) + NoLegend()
dev.off()
pdf('umap-features.pdf', width=8, height=8)
FeaturePlot(s_o, features=c("SoxN", "D", "pros"), ncol=2)
dev.off()
saveRDS(s_o, 's_o.bythebook.RDS')
