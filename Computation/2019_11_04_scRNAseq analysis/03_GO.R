library(dplyr)
library(tibble)
library(topGO)
library(org.Dm.eg.db)
library(rvest)
library(ggplot2)
setwd('/home/edsouza/Documents/Research/RussellLab/2019_11_04_scRNAseq analysis')


tf_list <- "https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/DNA_binding_proven_TF_maybe.html" %>%
  read_html() %>%
  html_nodes(xpath="/html/body/font/table") %>%
  html_table() %>% .[[1]] %>% dplyr::select(FBid) %>%
  mutate(symbol = mapIds(org.Dm.eg.db, 
                         keys=FBid, 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first"))


so_markers <- read.csv('../2020_04_29_ClusterMap/avalos.markers.csv')
fullGeneList <- readRDS('../2020_04_01_Avalos_scRNA/fullGeneList.RDS')



#library(ALL)
#data(ALL)
#data(geneList)

tblToNumeric <- function(tbl, bool=F) {
  nms <- names(tbl)
  if (bool) {
    out <- as.logical(tbl)
  } 
  else {
    out <- as.numeric(tbl)
  }
  names(out) <- nms
  out
}

geneSelFn <- function(clst, markerList=top100, manualList=NA) {
  if (is.na(manualList)) {
    genelist <- markerList %>% filter(cluster==clst) %>% .$gene %>% as.character()
  }
  else {
    genelist <- manualList
  }
  genes <- table(fullGeneList) %>% unlist()
  genes[names(genes) %in% genelist] <- T
  genes[!(names(genes) %in% genelist)] <- F
  tblToNumeric(genes, bool=F)
}
# selection <- function(value) {
#   value == 1
# }

selection <- function(bool) {return(bool)}


doAnalysis <- function(cluster, manualList=NA, manualListTitle=NA, markerList=top100) {
  
  # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
  bp <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Dm.eg.db", ID="symbol")
  GOdata <- new("topGOdata",
                ontology = "BP", 
                allGenes = geneSelFn(cluster, markerList, manualList),
                geneSel = selection,
                annot=annFUN.org,
                mapping = 'org.Dm.eg.db',
                ID = 'symbol')
  
  # TODO: see what else is always up or down when SoxN is up/down.
  # Look for specific gene names and find 
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  weightFisher <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
  #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  # https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
  allGO=usedGO(GOdata)
  all_res=GenTable(GOdata, classicFisher = resultFisher, 
                   weightFisher=weightFisher, 
                   orderBy='weightFisher', topNodes=length(allGO)) %>%
    mutate(padj = round(p.adjust(weightFisher, method="BH"),digits = 8))
  
  if (is.na(manualList)) {
    p1title <- paste('Cluster', cluster)
  } else {
    p1title <- paste('Manual list', manualListTitle)
  }
  p1 <- ggplot(all_res %>% filter(padj < 0.1), 
               aes(reorder(Term, -padj), 
                   -log10(padj),
                   fill=-log10(padj))) + 
    geom_bar(stat='identity') + 
    geom_hline(yintercept=-log10(0.1), linetype='dashed') +
    geom_hline(yintercept=-log10(0.05), linetype='dotted') +
    coord_flip() +
    labs(y='-log10(BH padj)', x='GO category', title=p1title) +
    theme(legend.position = "none")
  
  
  # allRes <- GenTable(GOdata, classicFisher = resultFisher,
  #                    weightFisher = weightFisher,
  #                    classicKS = resultKS, elimKS = resultKS.elim,
  #                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
  # allRes
  
  pvalFis <- score(weightFisher)
  pvalWeight <- score(weightFisher, whichGO = names(pvalFis))
  cor(pvalFis, pvalWeight)
  geneData(weightFisher)
  print(p1)
  
  out <- list(all_res)
  return(out)
}


# From ClusterMap, weak 0.04 assoc. with Avalos dataset clusters 6-8
doAnalysis(10, markerList = so_markers)
doAnalysis(0, markerList = so_markers,
           manualList = so_markers %>% 
             filter(cluster==10) %>% filter(gene %in% tf_list$symbol) 
           %>% .$gene,
           manualListTitle = 'Cluster 10 TFs')
# E(spl)mgamma-HLH	E(spl)m8-HLH	E(spl)m3-HLH	aop	ovo	pnt	**SoxN**	mod	Blimp-1	dmrt99B	
# Dek	E(spl)m5-HLH	bowl	E(spl)mdelta-HLH	E(spl)m7-HLH	tll	ci	ash2	CG9650
