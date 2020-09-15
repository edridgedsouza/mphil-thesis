library(dplyr)
library(ggplot2)


normalize <- T
if (normalize) {
  cormatfile <- 'cormat_lognorm.RDS'; distmatfile <- 'distmat_lognorm.RDS'
} else {
  cormatfile <- 'cormat.RDS'; distmatfile <- 'distmat.RDS'
}



# Load and get gene lists
cormat <- readRDS(cormatfile); distmat <- readRDS(distmatfile) %>% as.matrix()
Sox_closest <- distmat[,'SoxN'] %>% sort() %>% head(500) # 500 lowest Euclidean distance genes
Sox_correlated <- cormat[,'SoxN'] %>% sort() %>% tail(500) # 500 most correlated genes to SoxN

D_closest <- distmat[,'D'] %>% sort() %>% head(500)
D_correlated <- cormat[,'D'] %>% sort() %>% tail(500)

Sox_related <- union(names(Sox_closest), names(Sox_correlated))
D_related <- union(names(D_closest), names(D_correlated))
# Now plot again with actual cell type names
# new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(zinzen)
# zinzen <- RenameIdents(zinzen, new.cluster.ids)
# DimPlot(zinzen, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


qplot(names(Sox_correlated %>% head(20)), 
      Sox_correlated %>% head(20), 
      xlab='Gene', ylab='Correlation', main='Sox Correlated') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


########################
library(org.Dm.eg.db)
library(clusterProfiler)
# library(InterMineR)


genes.sox <- bitr(Sox_related, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Dm.eg.db, drop=F)
genes.d <- bitr(D_related, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Dm.eg.db, drop=F)


ego.sox <- enrichGO(gene      = genes.sox$ENTREZID,
                    ont           = "BP",
                    OrgDb = org.Dm.eg.db,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego.d <- enrichGO(gene      = genes.d$ENTREZID,
                  ont           = "BP",
                  OrgDb = org.Dm.eg.db,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

dotplot(ego.sox, showCategory=30)
dotplot(ego.d, showCategory=30)

cnetplot(ego.sox)
cnetplot(ego.d)

plotGOgraph(ego.sox)
plotGOgraph(ego.d)


library(UpSetR)
UpSetR::upset( fromList(input = list(D_closest = names(D_closest),
                                     D_correlated = names(D_correlated),
                                     Sox_closest = names(Sox_closest),
                                     Sox_correlated = names(Sox_correlated)
)
)
)
