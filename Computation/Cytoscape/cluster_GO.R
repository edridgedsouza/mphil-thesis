library(dplyr)
library(tibble)
library(topGO)
library(org.Dm.eg.db)
library(rvest)
library(ggplot2)
library(UpSetR)
library(purrr)

setwd('/home/edsouza/Documents/Research/RussellLab/Cytoscape/')


clusterLists <- readxl::read_excel('./SoxN D Two Degrees of Separation Nodes.csv')

so_markers <- read.table('../2020_07_03_PutItTogether/allen/so_markers.tsv', sep='\t')
fullGeneList <- readRDS('../2020_07_03_PutItTogether/fullGeneList.RDS')
tf_list <- "https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/DNA_binding_proven_TF_maybe.html" %>%
  read_html() %>%
  html_nodes(xpath="/html/body/font/table") %>%
  html_table() %>% .[[1]] %>% dplyr::select(FBid, 'Symbol/Name') %>%
  mutate(symbol = mapIds(org.Dm.eg.db, 
                         keys=FBid, 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first"))
tf_list$symbol[is.na(tf_list$symbol)] <- tf_list$`Symbol/Name`[is.na(tf_list$symbol)] %>% gsub('(.*) /.*','\\1', .)
tf_list$symbol[tf_list$symbol=='/'] <- NA # When org.Dm.eg.db can't convert IDs, use the table's annotation instead
SoxN_bind <- readxl::read_excel('../2020_07_03_PutItTogether/2020_04_16_SoxN_Dichaete.xlsx', sheet=1, col_names=T)
SoxN_target <- readxl::read_excel('../2020_07_03_PutItTogether/2020_04_16_SoxN_Dichaete.xlsx', sheet=2, col_names=T)
D_bind <- readxl::read_excel('../2020_07_03_PutItTogether/2020_04_16_SoxN_Dichaete.xlsx', sheet=3, col_names=T)
D_target <- readxl::read_excel('../2020_07_03_PutItTogether/2020_04_16_SoxN_Dichaete.xlsx', sheet=4, col_names=T)


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

geneSelFn <- function(clst, markerList=so_markers, manualList=NA) {
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


selection <- function(bool) {return(bool)}


doAnalysis <- function(cluster, manualList=NA, manualListTitle=NA, markerList=so_markers) {
  
  # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
  bp <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Dm.eg.db", ID="symbol")
  
  selectedGenes <- geneSelFn(cluster, markerList, manualList)
  numberOfGenes <- length(selectedGenes[selectedGenes != 0])
  GOdata <- new("topGOdata",
                ontology = "BP", 
                allGenes = selectedGenes,
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
  p1title <- paste0(p1title, ' (n = ', numberOfGenes, ')')
  
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
  #return(out)
  print(out)
  return(p1)
}
# Sometimes these intersections don't have enough genes and will crash
doAnalysis2 <- possibly(doAnalysis, otherwise=print('Not enough genes for comparison'))

savefig <- function(plotobj, out) {
  Cairo::CairoPDF(out, width=18, height=12)
  print(plotobj)
  dev.off()
}


doAnalysis2(cluster=0, manualList=clusterLists$`SoxN-D Two Degrees of Separation`,
           manualListTitle = 'SoxN-D Two Degrees of Separation') %>%
  savefig('figs/SoxN D Two Degrees of Separation.pdf')

doAnalysis2(cluster=0, manualList=clusterLists$`Top Cluster`,
            manualListTitle = 'Top Cluster') %>%
  savefig('figs/Top Cluster.pdf')

doAnalysis2(cluster=0, manualList=clusterLists$`Right Cluster`,
            manualListTitle = 'Right Cluster') %>%
  savefig('figs/Right Cluster.pdf')

doAnalysis2(cluster=0, manualList=clusterLists$`Bottom Cluster`,
            manualListTitle = 'Bottom Cluster') %>%
  savefig('figs/Bottom Cluster.pdf')
