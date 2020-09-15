library(dplyr)
library(tibble)
library(Seurat)

setwd('~/Documents/Research/RussellLab/2019_11_04_scRNAseq analysis/')
tables <- list()
for (i in 3:7) {
  datafile <- Sys.glob(paste0('data/*rep',i,'*'))
  tables[[i]] <- read.table(gzfile(datafile), sep='\t', header=T, quote="") %>% # Someon made a gene name with an apostrophe in it
    column_to_rownames('GENE') %>%
    CreateSeuratObject(project=paste0('rep',i))
}
# raw <- read.table(...., header=T) %>% column_to_rownames('GENE') 
zinzen <- merge(tables[[3]], y=c(tables[[4]], tables[[5]], tables[[6]], tables[[7]]),
                add.cell.ids= paste0('rep', 3:7),
                project="Zinzen")
rm(list='tables')

#################################### https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# nFeature_RNA: the number of genes
# nCount_RNA: the number of reads

replaceNA <- function(vec) {
  x <- vec
  x[is.na(x)] <- 0
  x
}

# For each cell, what % of reads map to genes fitting this pattern?
zinzen[["percent.trna"]] <- PercentageFeatureSet(zinzen, pattern = "^tRNA") %>% replaceNA
zinzen[["percent.rrna"]] <- PercentageFeatureSet(zinzen, pattern = "rRNA") %>% replaceNA
zinzen[["percent.mitoch"]] <- PercentageFeatureSet(zinzen, pattern = "^mt") %>% replaceNA
zinzen[["percent.mirna"]] <- PercentageFeatureSet(zinzen, pattern = "^mir") %>% replaceNA
zinzen[["percent.snrna"]] <- PercentageFeatureSet(zinzen, pattern = "snRNA") %>% replaceNA
zinzen[["percent.snorna"]] <- PercentageFeatureSet(zinzen, pattern = "snoRNA") %>% replaceNA
VlnPlot(zinzen, features = c("nFeature_RNA", "nCount_RNA", "percent.mitoch", "percent.trna", 
                             "percent.mirna",  "percent.rrna", "percent.snrna", "percent.snorna"), pt.size = 0.1)


plot1 <- FeatureScatter(zinzen, feature1 = "nCount_RNA", feature2 = "percent.trna")
plot2 <- FeatureScatter(zinzen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(zinzen, feature1 = "nCount_RNA", feature2 = "percent.rrna")
plot4 <- FeatureScatter(zinzen, feature1 = "percent.snrna", feature2 = "percent.snorna")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))


zinzen <- NormalizeData(zinzen, normalization.method = "LogNormalize")
zinzen <- FindVariableFeatures(zinzen, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(zinzen), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(zinzen)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


all.genes <- rownames(zinzen)
zinzen <- ScaleData(zinzen, features = all.genes)
zinzen <- RunPCA(zinzen, features = VariableFeatures(object = zinzen))
# print(zinzen[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(zinzen, dims = 6, reduction = "pca") # SoxN is top for dim 5, Sox21a to dim 6
# DimPlot(zinzen, reduction = "pca")



# Find dimensionality of dataset with heatmap inspection, JackStraw resampling, or elbow plot
# DimHeatmap(zinzen, dims = 1:6, cells = 1000, balanced = TRUE)
# zinzen <- JackStraw(zinzen, num.replicate = 100)
# zinzen <- ScoreJackStraw(zinzen, dims = 1:20)
# JackStrawPlot(zinzen, dims = 1:20)
# ElbowPlot(zinzen)
# Seems like first 20 dimensions are all meaningful

# Now cluster
zinzen <- FindNeighbors(zinzen, dims = 1:20) # Change dimensionality to most significant
zinzen <- FindClusters(zinzen, resolution = 0.5) # Adjust resolution parameter
# Find cluster values with Idents(zinzen)

# UMAP clustering
zinzen <- RunUMAP(zinzen, dims = 1:20)
Cairo::CairoPDF('figs/umap.pdf', width=10, height=10)
DimPlot(zinzen, reduction = "umap", label=T)
dev.off()
# tSNE clustering
# zinzen <- RunTSNE(zinzen, dims.use = 1:20, max_iter=2000, perplexity=5, check_duplicates = FALSE)
# DimPlot(zinzen, reduction = "tsne", label=T)
saveRDS(object=zinzen, file='zinzen.RDS')

#######################
# Way to correlate clusters with de novo marker identification
# Use this to see which clusters express SoxN, D
zinzen.markers <- FindAllMarkers(zinzen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GroupMarkers <- zinzen.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


# SoxN, D, vnd
# Known neuro markers: Prospero (pros), wor
FeaturePlot(zinzen, features = c("SoxN", "D", "vnd", "wor", "pros"))
FeaturePlot(zinzen, features = c("wor"))
RidgePlot(zinzen, features = c("SoxN", "D", "d", "vnd", "wor", "pros"), ncol = 2)
RidgePlot(zinzen, features = c("SoxN", "D",  "vnd", "pros"), ncol = 2)


top10 <- zinzen.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(zinzen, features = top10$gene) + NoLegend() # Use prior knowledge to assign cell type labels to each cluster

write.table(zinzen.markers, 'so_markers.tsv', sep='\t')

####################################################################################


library(ggplot2)
RidgePlot(zinzen, features = c("SoxN", "D"), ncol = 2)
# This gives the un-normalized data
# expressionClusters <- data.frame(SoxN = zinzen@assays$RNA@counts['SoxN',],
#                                  D = zinzen@assays$RNA@counts['D',],
#                                  cluster = Idents(zinzen))
# This gives the log-normalized data, as shown in the RidgePlot
vars <- c('SoxN', 'D')
expressionClusters <- FetchData(zinzen, slot="data", vars=vars)[, vars, drop=F] %>%
  mutate(cluster = Idents(zinzen), cell=colnames(zinzen))

ggplot(expressionClusters, aes(x=cluster, y = SoxN, fill=cluster)) + 
  geom_violin()
ggplot(expressionClusters, aes(x=cluster, y = D, fill=cluster)) + 
  geom_violin()


# Distance metrics to see what's similar to SoxN, D
# First create a subset of genes with at least 50 total reads across all samples
# in order to make the math less intense

normalize <- T
if (normalize) {
  subset <- FetchData(zinzen, slot='data', vars=rownames(zinzen)) %>% t() %>%  as.data.frame() %>%
    mutate(sums = rowSums(.), gene=rownames(.)) %>%
    filter(sums > 15) %>%
    select(-sums) %>% 
    column_to_rownames('gene')
  geneNames <- rownames(subset)
  subset <- subset %>% as.matrix()
  rownames(subset) <- geneNames
  names(subset) <- colnames(subset)
  
  distmat <- parallelDist::parDist(subset, threads=8) %>% as.matrix()
  saveRDS(distmat, 'distmat_lognorm.RDS'); rm(list='distmat')
  cormat <- cor(subset %>% t)
  saveRDS(cormat, 'cormat_lognorm.RDS'); rm(list='cormat')
  beepr::beep(8)
} else {
  subset <- zinzen@assays$RNA@counts %>% as.data.frame() %>%
    mutate(sums = rowSums(.), gene=rownames(.)) %>%
    filter(sums > 50) %>%
    select(-sums) %>% 
    column_to_rownames('gene')
  geneNames <- rownames(subset)
  subset <- subset %>% as.matrix()
  rownames(subset) <- geneNames
  names(subset) <- colnames(subset)
  
  distmat <- parallelDist::parDist(subset, threads=8) %>% as.matrix()
  saveRDS(distmat, 'distmat.RDS'); rm(list='distmat')
  cormat <- cor(subset %>% t)
  saveRDS(cormat, 'cormat.RDS'); rm(list='cormat')
  beepr::beep(8)
}

