library(dplyr)
library(tibble)
library(Seurat)

setwd('/home/edsouza/Documents/Research/RussellLab/2020_04_29_Allen_VentralNerveCord/')


# Reading allen dataset
type <- 'male' # 'male' or 'female
tables <- list()
for (i in 1:2) {
  datafile <- Sys.glob(paste0('data/*_', type, '*rep',i,'*'))
  tables[[i]] <- read.table(gzfile(datafile), sep='\t', header=T, quote="") %>% # Someone made a gene name with an apostrophe in it
    column_to_rownames('GENE') %>%
    CreateSeuratObject(project=paste0('rep',i))
}
s_o <- merge(tables[[1]], y=c(tables[[2]]),
                   add.cell.ids= paste0('rep', 1:2),
                   project=paste("Allen", type))
rm(list='tables'); gc()
saveRDS(s_o, 's_o.RDS')
# s_o <- readRDS('s_o.RDS')


#################################### https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# nFeature_RNA: the number of genes
# nCount_RNA: the number of reads

replaceNA <- function(vec) {
  x <- vec
  x[is.na(x)] <- 0
  x
}

# For each cell, what % of reads map to genes fitting this pattern?
s_o[["percent.trna"]] <- PercentageFeatureSet(s_o, pattern = "^tRNA") %>% replaceNA
s_o[["percent.rrna"]] <- PercentageFeatureSet(s_o, pattern = "rRNA") %>% replaceNA
s_o[["percent.mitoch"]] <- PercentageFeatureSet(s_o, pattern = "^mt") %>% replaceNA
s_o[["percent.mirna"]] <- PercentageFeatureSet(s_o, pattern = "^mir") %>% replaceNA
s_o[["percent.snrna"]] <- PercentageFeatureSet(s_o, pattern = "snRNA") %>% replaceNA
s_o[["percent.snorna"]] <- PercentageFeatureSet(s_o, pattern = "snoRNA") %>% replaceNA
VlnPlot(s_o, features = c("nFeature_RNA", "nCount_RNA", "percent.mitoch", "percent.trna", 
                          "percent.mirna",  "percent.rrna", "percent.snrna", "percent.snorna"), 
        ncol=4, pt.size = 0.1)


plot1 <- FeatureScatter(s_o, feature1 = "nCount_RNA", feature2 = "percent.trna")
plot2 <- FeatureScatter(s_o, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(s_o, feature1 = "nCount_RNA", feature2 = "percent.rrna")
plot4 <- FeatureScatter(s_o, feature1 = "percent.snrna", feature2 = "percent.snorna")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))


s_o <- NormalizeData(s_o, normalization.method = "LogNormalize")
s_o <- FindVariableFeatures(s_o, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(s_o), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(s_o)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


all.genes <- rownames(s_o)
s_o <- ScaleData(s_o, features = all.genes) 
  # saveRDS(s_o, 's_o.scaled.RDS')
s_o <- RunPCA(s_o, features = VariableFeatures(object = s_o))
# print(s_o[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(s_o, dims = 6, reduction = "pca") # SoxN is top for dim 5, Sox21a to dim 6
# DimPlot(s_o, reduction = "pca")



# Find dimensionality of dataset with heatmap inspection, JackStraw resampling, or elbow plot
# DimHeatmap(s_o, dims = 1:6, cells = 1000, balanced = TRUE)
# s_o <- JackStraw(s_o, num.replicate = 100)
# s_o <- ScoreJackStraw(s_o, dims = 1:20)
# JackStrawPlot(s_o, dims = 1:20)
# ElbowPlot(s_o)
# Seems like first 20 dimensions are all meaningful

# Now cluster
s_o <- FindNeighbors(s_o, dims = 1:20) # Change dimensionality to most significant
s_o <- FindClusters(s_o, resolution = 0.5) # Adjust resolution parameter
# Find cluster values with Idents(s_o)

# UMAP clustering
s_o <- RunUMAP(s_o, dims = 1:20)
DimPlot(s_o, reduction = "umap", label=T)

# tSNE clustering
# s_o <- RunTSNE(s_o, dims.use = 1:20, max_iter=2000, perplexity=5, check_duplicates = FALSE)
# DimPlot(s_o, reduction = "tsne", label=T)

#######################
# Way to correlate clusters with de novo marker identification
# Use this to see which clusters express SoxN, D
s_o.markers <- FindAllMarkers(s_o, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GroupMarkers <- s_o.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


# SoxN, D, vnd
# Known neuro markers: Prospero (pros), wor
FeaturePlot(s_o, features = c("SoxN", "D", "vnd", "wor", "pros"))
FeaturePlot(s_o, features = c("wor"))
RidgePlot(s_o, features = c("SoxN", "D", "d", "vnd", "wor", "pros"), ncol = 2)
RidgePlot(s_o, features = c("SoxN", "D",  "vnd", "pros"), ncol = 2)


top100 <- s_o.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.table(s_o.markers, 'so_markers.tsv', sep='\t')
write.table(top100, 'top100markers.tsv', sep='\t')
write.table(GroupMarkers, 'top10markers.tsv', sep='\t')
DoHeatmap(s_o, features = top10$gene) + NoLegend() # Use prior knowledge to assign cell type labels to each cluster


####################################################################################


library(ggplot2)
RidgePlot(s_o, features = c("SoxN", "D"), ncol = 2)
# This gives the un-normalized data
# expressionClusters <- data.frame(SoxN = s_o@assays$RNA@counts['SoxN',],
#                                  D = s_o@assays$RNA@counts['D',],
#                                  cluster = Idents(s_o))
# This gives the log-normalized data, as shown in the RidgePlot
vars <- c('SoxN', 'D')
expressionClusters <- FetchData(s_o, slot="data", vars=vars)[, vars, drop=F] %>%
  mutate(cluster = Idents(s_o), cell=colnames(s_o))

ggplot(expressionClusters, aes(x=cluster, y = SoxN, fill=cluster)) + 
  geom_violin()
ggplot(expressionClusters, aes(x=cluster, y = D, fill=cluster)) + 
  geom_violin()


# Distance metrics to see what's similar to SoxN, D
# First create a subset of genes with at least 50 total reads across all samples
# in order to make the math less intense

normalize <- T
if (normalize) {
  subset <- FetchData(s_o, slot='data', vars=rownames(s_o)) %>% t() %>%  as.data.frame() %>%
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
  subset <- s_o@assays$RNA@counts %>% as.data.frame() %>%
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

