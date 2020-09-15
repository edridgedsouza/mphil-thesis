library(dplyr)
library(tibble)
library(topGO)
library(org.Dm.eg.db)
library(rvest)
library(ggplot2)
setwd('/home/edsouza/Documents/Research/RussellLab/2020_04_01_Avalos_scRNA')


tf_list <- "https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/DNA_binding_proven_TF_maybe.html" %>%
  read_html() %>%
  html_nodes(xpath="/html/body/font/table") %>%
  html_table() %>% .[[1]] %>% dplyr::select(FBid) %>%
  mutate(symbol = mapIds(org.Dm.eg.db, 
                keys=FBid, 
                column="SYMBOL", 
                keytype="ENSEMBL",
                multiVals="first"))


so_markers <- read.table('so_markers.tsv', sep='\t')
top100 <- read.table('top100markers.tsv', sep='\t')
top10 <- read.table('top10markers.tsv', sep='\t')
fullGeneList <- readRDS('fullGeneList.RDS')



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


doAnalysis(6, markerList = so_markers)
doAnalysis(0, markerList = so_markers,
           manualList = so_markers %>% 
             filter(cluster==6) %>% filter(gene %in% tf_list$symbol) 
           %>% .$gene,
           manualListTitle = 'Cluster 6 TFs')
doAnalysis(7, markerList = so_markers)
doAnalysis(0, markerList = so_markers,
           manualList = so_markers %>% 
             filter(cluster==7) %>% 
             filter(gene %in% tf_list$symbol) %>% .$gene,
           manualListTitle = 'Cluster 7 TFs')
doAnalysis(8, markerList = so_markers)
doAnalysis(0, markerList = so_markers,
           manualList = so_markers %>% 
             filter(cluster==8) %>% 
             filter(gene %in% tf_list$symbol) %>% .$gene,
           manualListTitle = 'Cluster 8 TFs')

doAnalysis(13, markerList = so_markers)
doAnalysis(0, markerList = so_markers,
           manualList = so_markers %>% 
             filter(cluster==13) %>% 
             filter(gene %in% tf_list$symbol) %>% .$gene,
           manualListTitle = 'Cluster 13 TFs')

doAnalysis(0, manualList = so_markers %>% filter(cluster %in% c(6,7,8)) %>% .$gene %>% as.character(),
           manualListTitle = 'Cluster 6-8 all')
doAnalysis(0, manualList = so_markers %>% 
             filter(cluster %in% c(6,7,8)) %>% 
             filter(gene %in% tf_list$symbol) %>% .$gene,
           manualListTitle = 'Cluster 6-8 TFs')

top100 %>% 
  filter(cluster %in% c(6,7,8)) %>% 
  filter(gene %in% tf_list$symbol)  -> cluster678_tfs


doAnalysis(0, manualList = so_markers %>% filter(cluster %in% c(6,7,8, 10)) %>% .$gene %>% as.character(),
           manualListTitle = 'Cluster 6-8, 10 all')
doAnalysis(0, manualList = so_markers %>% 
             filter(cluster %in% c(6,7,8, 10)) %>% 
             filter(gene %in% tf_list$symbol) %>% .$gene,
           manualListTitle = 'Cluster 6-8, 10 TFs')
top100 %>% 
  filter(cluster %in% c(6,7,8, 10)) %>% 
  filter(gene %in% tf_list$symbol)  -> cluster678_10_tfs


library(UpSetR)
list(#tfs = tf_list$symbol,
     cluster6 = so_markers %>% filter(cluster==6) %>% .$gene,
     cluster7 = so_markers %>% filter(cluster==7) %>% .$gene,
     cluster8 = so_markers %>% filter(cluster==8) %>% .$gene,
     cluster10 = so_markers %>% filter(cluster==10) %>% .$gene) %>%
  fromList() %>%
  upset(empty.intersections = "on", order.by = "freq", nintersects=100)

# Common factors in these 4 categories
so_markers %>% filter(cluster==6) %>% .$gene %>%
  intersect(so_markers %>% filter(cluster==7) %>% .$gene) %>%
  intersect(so_markers %>% filter(cluster==8) %>% .$gene) %>%
  intersect(so_markers %>% filter(cluster==10) %>% .$gene) -> SoxN_assoc
doAnalysis(0, manualList = SoxN_assoc,
           manualListTitle = 'SoxN Associated Cluster Intersection Consensus')
SoxN_assoc[SoxN_assoc %in% tf_list$symbol] # Intersection TF consensus, too small to do GO on
top100 %>% 
  filter(cluster %in% c(6,7,8)) %>% .$gene %>%
  doAnalysis(0, manualList = ., manualListTitle='SoxN Associated Cluster Union Consensus')
cluster678_10_tfs$gene %>%
  doAnalysis(0, manualList = ., manualListTitle='SoxN Associated Cluster Union TF Consensus')
# fruitless (fru) encodes a BTB zinc finger transcription factor that contributes to sexual differentiation of the neural circuits underlying male sexual behavior.
# castor (cas) encodes a transcription factor expressed in the latest stage of embryonic neuroblast lineages
# modulo (mod) encodes a nucleolar protein that is required for growth of proliferative cells through its association with the proto-oncogene Myc.
# klumpfuss (klu) zinc finger - early growth response family - involved in determination of the identities of secondary precursor cells within neuroblast lineages in the central nervous system 

# Compare sets with validated list of SoxN/D binding targets / regulatory targets
SoxN_bind <- readxl::read_excel('suppData/2020_04_16_SoxN_Dichaete.xlsx', sheet=1, col_names=T)
SoxN_target <- readxl::read_excel('suppData/2020_04_16_SoxN_Dichaete.xlsx', sheet=2, col_names=T)

list(#tfs = tf_list$symbol,
  cluster6 = so_markers %>% filter(cluster==6) %>% .$gene,
  cluster7 = so_markers %>% filter(cluster==7) %>% .$gene,
  cluster8 = so_markers %>% filter(cluster==8) %>% .$gene,
  cluster10 = so_markers %>% filter(cluster==10) %>% .$gene, 
  SoxN_bind = SoxN_bind$Symbol) %>% # SoxN_bind$Symbol
  fromList() %>%
  upset(empty.intersections = "on", order.by = "freq", nintersects=100)

so_markers %>% filter(cluster==6) %>% .$gene %>%
  intersect(so_markers %>% filter(cluster==7) %>% .$gene) %>%
  intersect(so_markers %>% filter(cluster==8) %>% .$gene) %>%
  intersect(so_markers %>% filter(cluster==10) %>% .$gene) %>%
  intersect(SoxN_target$Symbol)

so_markers %>% filter(cluster %in% c(6,7,8, 10)) %>% .$gene %>%
  intersect(SoxN_target$Symbol) %>%
  doAnalysis(0, manualList = ., manualListTitle='Clusters union 6-8,10 targeted by SoxN')
# "CG17265" "klu"     "cas"     "Chd64"   "RpS15"
# "CG17265" is Ccdc85, involved with RNAi.
#  Chd64 (Chd64) encodes a protein involved in the juvenile hormone mediated signaling pathway.
# Ribosomal protein S15 (RpS15) encodes a protein that plays a role in assembly, structure and function of the ribosome and functions in protein synthesis.


top100 %>% mutate(TF = gene %in% tf_list$symbol) %>%
  write.table('top100markers.tsv', sep='\t')
so_markers %>% mutate(TF = gene %in% tf_list$symbol) %>%
  write.table('allmarkers.tsv', sep='\t')



# _____________________________________________________________________________________
so_markers %>% filter(cluster %in% c(6,7,8,10,13)) %>% .$gene %>%
  intersect(SoxN_bind$Symbol) %>%
  cat(sep='\n') # Put into Cytoscape


so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene %>%
  intersect(SoxN_bind$Symbol) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(7)) %>% .$gene %>%
  intersect(SoxN_bind$Symbol) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(7)) %>% .$gene %>%
  setdiff(so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(7)) %>% .$gene %>%
  setdiff(so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene) %>%
  intersect(SoxN_bind$Symbol) %>%
  cat(sep='\n')

# so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene %>%
#   setdiff(so_markers %>% filter(cluster %in% c(7)) %>% .$gene) %>%
#   cat(sep='\n')

so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene %>%
  setdiff(so_markers %>% filter(cluster %in% c(7)) %>% .$gene) %>%
  intersect(SoxN_bind$Symbol) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene %>%
  intersect(SoxN_target$Symbol) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(7)) %>% .$gene %>%
  intersect(SoxN_target$Symbol) %>%
  cat(sep='\n')

so_markers %>% filter(cluster %in% c(7)) %>% .$gene %>%
  setdiff(so_markers %>% filter(cluster %in% c(6,8)) %>% .$gene) %>%
  intersect(SoxN_target$Symbol) %>%
  cat(sep='\n')
