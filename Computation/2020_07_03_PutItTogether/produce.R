library(dplyr)
library(tibble)
library(topGO)
library(org.Dm.eg.db)
library(rvest)
library(ggplot2)
library(UpSetR)
library(purrr)
setwd('/home/edsouza/Documents/Research/RussellLab/2020_07_03_PutItTogether/')

# Starting setup: 3 folders named zinzen, avalos, allen. Each have so_markers.tsv inside
# which was obtained from the original Seurat analysis.
# Also included in main dir: 2020_04_16_SoxN_Dichaete.xlsx, fullGeneList.RDS, TableS1.xlsx

tf_list <-
  "https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/DNA_binding_proven_TF_maybe.html" %>%
  read_html() %>%
  html_nodes(xpath = "/html/body/font/table") %>%
  html_table() %>% .[[1]] %>% dplyr::select(FBid, 'Symbol/Name') %>%
  mutate(
    symbol = mapIds(
      org.Dm.eg.db,
      keys = FBid,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
tf_list$symbol[is.na(tf_list$symbol)] <-
  tf_list$`Symbol/Name`[is.na(tf_list$symbol)] %>% gsub('(.*) /.*', '\\1', .)
tf_list$symbol[tf_list$symbol == '/'] <-
  NA # When org.Dm.eg.db can't convert IDs, use the table's annotation instead


SoxN_bind <-
  readxl::read_excel('2020_04_16_SoxN_Dichaete.xlsx',
                     sheet = 1,
                     col_names = T)
SoxN_target <-
  readxl::read_excel('2020_04_16_SoxN_Dichaete.xlsx',
                     sheet = 2,
                     col_names = T)
D_bind <-
  readxl::read_excel('2020_04_16_SoxN_Dichaete.xlsx',
                     sheet = 3,
                     col_names = T)
D_target <-
  readxl::read_excel('2020_04_16_SoxN_Dichaete.xlsx',
                     sheet = 4,
                     col_names = T)
fullGeneList <- readRDS('fullGeneList.RDS')

createMarkerFile <- function(dir) {
  so_markers <- read.table(paste0(dir, '/so_markers.tsv'), sep = '\t')
  so_markers <- so_markers %>% mutate(
    TF = gene %in% tf_list$symbol,
    SoxN_bind = gene %in% SoxN_bind$Symbol,
    SoxN_target = gene %in% SoxN_target$Symbol,
    D_bind = gene %in% D_bind$Symbol,
    D_target = gene %in% D_target$Symbol
  )
  
  write.table(
    so_markers,
    paste0(dir, '/', dir, '_markers.tsv'),
    sep = '\t',
    row.names = F
  )
}

savefig <- function(plotobj,
                    out,
                    width = 7.5,
                    height = 7.5) {
  Cairo::CairoPDF(out, width = width, height = height)
  print(plotobj)
  dev.off()
}

createPlots <- function(dir) {
  so_markers <- read.table(paste0(dir, '/so_markers.tsv'), sep = '\t')
  hotClusters <-
    so_markers[so_markers$gene %in% c('SoxN', 'D'), ]$cluster
  soxnClusters <-
    so_markers[so_markers$gene %in% c('SoxN'), ]$cluster
  dClusters <- so_markers[so_markers$gene %in% c('D'), ]$cluster
  
  
  
  tblToNumeric <- function(tbl, bool = F) {
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
  
  geneSelFn <- function(clst,  markerList = so_markers, manualList = NA) {
      if (is.na(manualList)) {
        genelist <-
          markerList %>% filter(cluster == clst) %>% .$gene %>% as.character()
      }
      else {
        genelist <- manualList
      }
      genes <- table(fullGeneList) %>% unlist()
      genes[names(genes) %in% genelist] <- T
      genes[!(names(genes) %in% genelist)] <- F
      tblToNumeric(genes, bool = F)
    }
  
  selection <- function(bool) {
    return(bool)
  }
  
  doAnalysis <-
    function(cluster, manualList = NA,  manualListTitle = NA, markerList = so_markers) {
      # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
      bp <-
        annFUN.org(
          whichOnto = "BP",
          feasibleGenes = NULL,
          mapping = "org.Dm.eg.db",
          ID = "symbol"
        )
      
      selectedGenes <- geneSelFn(cluster, markerList, manualList)
      numberOfGenes <- length(selectedGenes[selectedGenes != 0])
      GOdata <- new(
        "topGOdata",
        ontology = "BP",
        allGenes = selectedGenes,
        geneSel = selection,
        annot = annFUN.org,
        mapping = 'org.Dm.eg.db',
        ID = 'symbol'
      )
      
      # TODO: see what else is always up or down when SoxN is up/down.
      # Look for specific gene names and find
      resultFisher <-
        runTest(GOdata, algorithm = "classic", statistic = "fisher")
      weightFisher <-
        runTest(GOdata, algorithm = 'weight01', statistic = 'fisher')
      #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
      #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
      
      # https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
      allGO = usedGO(GOdata)
      all_res = GenTable(
        GOdata,
        classicFisher = resultFisher,
        weightFisher = weightFisher,
        orderBy = 'weightFisher',
        topNodes = length(allGO)
      ) %>%
        mutate(padj = round(p.adjust(weightFisher, method = "BH"), digits = 8))
      
      if (is.na(manualList)) {
        p1title <- paste('Cluster', cluster)
      } else {
        p1title <- paste('Manual list', manualListTitle)
      }
      p1title <- paste0(p1title, ' (n = ', numberOfGenes, ')')
      
      
      # Limit output to only 50 most significant if there are too many
      # Comment this block   out before putting full data dump on GitHub
      if (dim(all_res)[1] > 50) {
        all_res <- all_res %>% top_n(-50, padj)
      }
      
      
      p1 <- ggplot(all_res %>% filter(padj < 0.1),
                   aes(reorder(Term,-padj),-log10(padj),
                       fill = -log10(padj))) +
        geom_bar(stat = 'identity') +
        geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), linetype = 'dotted') +
        coord_flip() +
        labs(y = '-log10(BH padj)', x = 'GO category', title = p1title) +
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
      # return(out)
      return(p1)
    }
  # Sometimes these intersections don't have enough genes and will crash
  doAnalysis2 <- possibly(doAnalysis, otherwise = print('Not enough genes for comparison'))
  
  analyzeCluster <- function(clusterNumber) {
    doAnalysis2(clusterNumber, markerList = so_markers) %>%
      savefig(paste0(dir, '/cluster', clusterNumber, '.pdf'))
    doAnalysis2(
      0,
      markerList = so_markers,
      manualList = so_markers %>%
        filter(cluster == clusterNumber) %>% filter(gene %in% tf_list$symbol)
      %>% .$gene,
      manualListTitle = paste0('Cluster ', clusterNumber, ' TFs')
    ) %>%
      savefig(paste0(dir, '/cluster', clusterNumber, '_TFs.pdf'))
  }
  
  # Gives plots for union of all input clusters
  analyzeAllClusters <-  function(clusterNumberList, prefix = 'Combined') {
      if (length(clusterNumberList) == 0) {
        print('No elements in this list')
        return(0)
      }
      
      doAnalysis2(
        0,
        manualList = so_markers %>% filter(cluster %in% clusterNumberList) %>% .$gene %>% as.character(),
        manualListTitle = paste0(prefix, ' Cluster superset all')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterSuperset_all.pdf'))
      doAnalysis2(
        0,
        manualList = so_markers %>%
          filter(cluster %in% clusterNumberList) %>%
          filter(gene %in% tf_list$symbol) %>% .$gene,
        manualListTitle = paste0(prefix, ' Cluster superset TFs')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterSuperset_TFs.pdf'))
      doAnalysis2(
        0,
        manualList = so_markers %>%
          filter(cluster %in% clusterNumberList) %>%
          filter(gene %in% SoxN_bind$Symbol) %>% .$gene,
        manualListTitle = paste0(prefix, ' Cluster superset, SoxN bound')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterSuperset_SoxNBind.pdf'))
      doAnalysis2(
        0,
        manualList = so_markers %>%
          filter(cluster %in% clusterNumberList) %>%
          filter(gene %in% D_bind$Symbol) %>% .$gene,
        manualListTitle = paste0(prefix, ' Cluster superset, D bound')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterSuperset_DBind.pdf'))
    }
  
  # Give plots for INTERSECTION of all input clusters
  clusterIntersect <- function(clusterNumberList, prefix = 'Combined') {
      if (length(clusterNumberList) == 0) {
        print('No elements in this list')
        return(0)
      }
      
      
      intersectionList <- list()
      for (i in seq_along(clusterNumberList)) {
        intersectionList[[i]] = so_markers %>% filter(cluster == clusterNumberList[i]) %>% .$gene
      }
      intersectionOfAll <- Reduce(intersect, intersectionList)
      
      if (length(intersectionOfAll) == 0) {
        print('There are elements in this list, but their intersection is empty')
        return(0)
      }
      
      
      
      doAnalysis2(
        0,
        manualList = intersectionOfAll,
        manualListTitle =  paste0(prefix, ' Cluster intersection all')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterIntersection_all.pdf'))
      
      doAnalysis2(
        0,
        manualList = intersectionOfAll %>% intersect(tf_list$symbol),
        manualListTitle =  paste0(prefix, ' Cluster intersection TFs')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterIntersection_tfs.pdf'))
      
      doAnalysis2(
        0,
        manualList = intersectionOfAll %>% intersect(SoxN_bind$Symbol),
        manualListTitle =  paste0(prefix, ' Cluster intersection, SoxN bound')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterIntersection_SoxNBind.pdf'))
      
      doAnalysis2(
        0,
        manualList = intersectionOfAll %>% intersect(D_bind$Symbol),
        manualListTitle =  paste0(prefix, ' Cluster intersection, D bound')
      ) %>%
        savefig(paste0(dir, '/', prefix, '_clusterIntersection_DBind.pdf'))
      
    }
  
  ### Actually running now
  for (cluster in hotClusters) {
    print(paste('Analyzing', dir, 'cluster', cluster))
    analyzeCluster(cluster)
  }
  print('Running analyzeAllClusters block')
  print(paste(dir, 'union Combined'))
  analyzeAllClusters(hotClusters, prefix = 'Combined')
  print(paste(dir, 'union SoxN'))
  analyzeAllClusters(soxnClusters, prefix = 'SoxN')
  print(paste(dir, 'union D'))
  analyzeAllClusters(dClusters, prefix = 'D')
  
  print('Running clusterIntersect block')
  print(paste(dir, 'intersect Combined'))
  clusterIntersect(hotClusters, prefix = 'Combined')
  print(paste(dir, 'intersect SoxN'))
  clusterIntersect(soxnClusters, prefix = 'SoxN')
  print(paste(dir, 'intersect D'))
  clusterIntersect(dClusters, prefix = 'D')
  
  print('Running UpSet block')
  clusterLst <- list()
  for (i in hotClusters) {
    clusterLst[[paste0('cluster', i)]] = so_markers %>% filter(cluster == i) %>% .$gene
  }
  clusterLst %>% fromList() %>%
    upset(
      empty.intersections = "on",
      order.by = "freq",
      nintersects = 100,
      nsets = length(hotClusters)
    ) %>%
    savefig(paste0(dir, '/upset.pdf'),
            width = 6,
            height = 4)
  
  
  #### Ran union of SoxN+D clusters. What else?
  # Intersection of all SoxN? Of all D? Of both? Bit harder to code
}





############################################################
createMarkerFile('zinzen')
createMarkerFile('avalos')
createMarkerFile('allen')

createPlots('zinzen')
createPlots('avalos')
createPlots('allen') # Allen markers must be "by the book". Cluster 108 has both SoxN and D
# All GO plots are Biological Process


# Note: the GO plots don't tell you which gene has defined the cluster
# In most cases it's SoxN, but you should check in the spreadsheet to be sure.



###############################################################################################
# Ad-hoc analyses


####################################### Ad-hoc analysis for Allen cluster 108
markers <-
  read.table('allen/allen_markers.tsv', sep = '\t', header = T) %>% filter(cluster ==
                                                                             108)
cluster108 <- readxl::read_excel('TableS1.xlsx', sheet = 1)
# Need to replace cluster108 with the synonyms that are in the original marker list
# Cytoscape converted some gene names to synonyms
# "AcrE"     "dor"      NA         "bl"       "pUf68"    "Inr-a"
# "snRNP70K" "CG11984"  "CG42575"  "Mpcp"     "Hsp60" "CG7967"  "vfl"
# Yeah it's repetitive. Sue me
cluster108[cluster108 == 'AcrE'] <- 'nAChRalpha2'
cluster108[cluster108 == 'dor'] <- 'DOR'
cluster108[cluster108 == 'bl'] <- 'HnRNP-K'
cluster108[cluster108 == 'pUf68'] <- 'hfp'
cluster108[cluster108 == 'Inr-a'] <- 'Pcf11'
cluster108[cluster108 == 'snRNP70K'] <- 'snRNP-U1-70K'
cluster108[cluster108 == 'CG11984'] <- 'Kcmf1'
cluster108[cluster108 == 'CG42575'] <- 'NaPi-III'
cluster108[cluster108 == 'Mpcp'] <- 'Mpcp2'
cluster108[cluster108 == 'Hsp60'] <- 'Hsp60A'
cluster108[cluster108 == 'CG7967'] <- 'Vta1'
cluster108[cluster108 == 'vfl'] <- 'zld'


summary <- function(subset) {
  out <- subset %>%
    summarize(
      n = n(),
      TF = sum(TF),
      SoxN_bind = sum(SoxN_bind),
      SoxN_target = sum(SoxN_target),
      D_bind = sum(D_bind),
      D_target = sum(D_target)
    )
  
  out2 <- out %>% `/`(.$n)
  
  print(out)
  print(out2)
}

markers %>% summary()
markers[markers$gene %in% cluster108$`SoxN-D Two Degrees of Separation`, ] %>% summary()
markers[markers$gene %in% cluster108$`Cluster A`, ] %>% summary()
markers[markers$gene %in% cluster108$`Cluster B`, ] %>% summary()
markers[markers$gene %in% cluster108$`Cluster C`, ] %>% summary()

markers[markers$gene %in% unique(c(
  cluster108$`Cluster A`,
  cluster108$`Cluster B`,
  cluster108$`Cluster C`
)), ]  %>%  summary()



################################## Ad-hoc correction of the UpSet plot dimensions

so_markers <- read.table('allen/so_markers.tsv', sep = '\t')
clusterLst <- list()
for (i in c(135, 167, 126, 108, 53, 159, 151, 140)) {
  clusterLst[[paste0('cluster', i)]] = so_markers %>% filter(cluster == i) %>% .$gene
}
clusterLst %>% fromList() %>%
  upset(
    empty.intersections = "on",
    order.by = "freq",
    nintersects = 30,
    nsets = 8
  ) %>%
  savefig('allen/upset-adhoc.pdf',
          width = 6 * 1,
          height = 4 * 1)


####################################################### Ad-hoc zinzen cluster 12 for figure
# Re-copying all the helper functions within createPlots().
# Slightly tweaked them to work outside of the main function environment
# u mad, compsci majors who expect clean code? :)

tblToNumeric <- function(tbl, bool = F) {
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

geneSelFn <- function(clst,  markerList = so_markers, manualList = NA) {
  if (is.na(manualList)) {
    genelist <-
      markerList %>% filter(cluster == clst) %>% .$gene %>% as.character()
  }
  else {
    genelist <- manualList
  }
  genes <- table(fullGeneList) %>% unlist()
  genes[names(genes) %in% genelist] <- T
  genes[!(names(genes) %in% genelist)] <- F
  tblToNumeric(genes, bool = F)
}

selection <- function(bool) {
  return(bool)
}

doAnalysis <-
  function(cluster, manualList = NA,  manualListTitle = NA, markerList = so_markers) {
    # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
    bp <-
      annFUN.org(
        whichOnto = "BP",
        feasibleGenes = NULL,
        mapping = "org.Dm.eg.db",
        ID = "symbol"
      )
    
    selectedGenes <- geneSelFn(cluster, markerList, manualList)
    numberOfGenes <- length(selectedGenes[selectedGenes != 0])
    GOdata <- new(
      "topGOdata",
      ontology = "BP",
      allGenes = selectedGenes,
      geneSel = selection,
      annot = annFUN.org,
      mapping = 'org.Dm.eg.db',
      ID = 'symbol'
    )
    
    # TODO: see what else is always up or down when SoxN is up/down.
    # Look for specific gene names and find
    resultFisher <-
      runTest(GOdata, algorithm = "classic", statistic = "fisher")
    weightFisher <-
      runTest(GOdata, algorithm = 'weight01', statistic = 'fisher')
    #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    # https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
    allGO = usedGO(GOdata)
    all_res = GenTable(
      GOdata,
      classicFisher = resultFisher,
      weightFisher = weightFisher,
      orderBy = 'weightFisher',
      topNodes = length(allGO)
    ) %>%
      mutate(padj = round(p.adjust(weightFisher, method = "BH"), digits = 8))
    
    if (is.na(manualList)) {
      p1title <- paste('Cluster', cluster)
    } else {
      p1title <- paste('Manual list', manualListTitle)
    }
    p1title <- paste0(p1title, ' (n = ', numberOfGenes, ')')
    
    
    # Limit output to only 50 most significant if there are too many
    # Comment this block   out before putting full data dump on GitHub
    if (dim(all_res)[1] > 50) {
      all_res <- all_res %>% top_n(-50, padj)
    }
    
    
    p1 <- ggplot(all_res %>% filter(padj < 0.1),
                 aes(reorder(Term,-padj),-log10(padj),
                     fill = -log10(padj))) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted') +
      coord_flip() +
      labs(y = '-log10(BH padj)', x = 'GO category', title = p1title) +
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
    # return(out)
    return(p1)
  }
# Sometimes these intersections don't have enough genes and will crash
doAnalysis2 <- possibly(doAnalysis, otherwise = print('Not enough genes for comparison'))


### PART THAT I MODIFIED
analyzeCluster <- function(clusterNumber, markerList = so_markers, height=4, width=4) {
  doAnalysis2(clusterNumber, markerList = so_markers) %>%
    savefig(paste0(dir, '/cluster', clusterNumber, '.pdf'), height=height, width=width)
  doAnalysis2(
    0,
    markerList = so_markers,
    manualList = so_markers %>%
      filter(cluster == clusterNumber) %>% filter(gene %in% tf_list$symbol)
    %>% .$gene,
    manualListTitle = paste0('Cluster ', clusterNumber, ' TFs')
  ) %>%
    savefig(paste0(dir, '/cluster', clusterNumber, '_TFs.pdf'), height=height, width=width)
}

# Gives plots for union of all input clusters
analyzeAllClusters <-  function(clusterNumberList, prefix = 'Combined') {
  if (length(clusterNumberList) == 0) {
    print('No elements in this list')
    return(0)
  }
  
  doAnalysis2(
    0,
    manualList = so_markers %>% filter(cluster %in% clusterNumberList) %>% .$gene %>% as.character(),
    manualListTitle = paste0(prefix, ' Cluster superset all')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterSuperset_all.pdf'))
  doAnalysis2(
    0,
    manualList = so_markers %>%
      filter(cluster %in% clusterNumberList) %>%
      filter(gene %in% tf_list$symbol) %>% .$gene,
    manualListTitle = paste0(prefix, ' Cluster superset TFs')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterSuperset_TFs.pdf'))
  doAnalysis2(
    0,
    manualList = so_markers %>%
      filter(cluster %in% clusterNumberList) %>%
      filter(gene %in% SoxN_bind$Symbol) %>% .$gene,
    manualListTitle = paste0(prefix, ' Cluster superset, SoxN bound')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterSuperset_SoxNBind.pdf'))
  doAnalysis2(
    0,
    manualList = so_markers %>%
      filter(cluster %in% clusterNumberList) %>%
      filter(gene %in% D_bind$Symbol) %>% .$gene,
    manualListTitle = paste0(prefix, ' Cluster superset, D bound')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterSuperset_DBind.pdf'))
}

# Give plots for INTERSECTION of all input clusters
clusterIntersect <- function(clusterNumberList, prefix = 'Combined') {
  if (length(clusterNumberList) == 0) {
    print('No elements in this list')
    return(0)
  }
  
  
  intersectionList <- list()
  for (i in seq_along(clusterNumberList)) {
    intersectionList[[i]] = so_markers %>% filter(cluster == clusterNumberList[i]) %>% .$gene
  }
  intersectionOfAll <- Reduce(intersect, intersectionList)
  
  if (length(intersectionOfAll) == 0) {
    print('There are elements in this list, but their intersection is empty')
    return(0)
  }
  
  
  
  doAnalysis2(
    0,
    manualList = intersectionOfAll,
    manualListTitle =  paste0(prefix, ' Cluster intersection all')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterIntersection_all.pdf'))
  
  doAnalysis2(
    0,
    manualList = intersectionOfAll %>% intersect(tf_list$symbol),
    manualListTitle =  paste0(prefix, ' Cluster intersection TFs')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterIntersection_tfs.pdf'))
  
  doAnalysis2(
    0,
    manualList = intersectionOfAll %>% intersect(SoxN_bind$Symbol),
    manualListTitle =  paste0(prefix, ' Cluster intersection, SoxN bound')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterIntersection_SoxNBind.pdf'))
  
  doAnalysis2(
    0,
    manualList = intersectionOfAll %>% intersect(D_bind$Symbol),
    manualListTitle =  paste0(prefix, ' Cluster intersection, D bound')
  ) %>%
    savefig(paste0(dir, '/', prefix, '_clusterIntersection_DBind.pdf'))
  
}

savefig <- function(plotobj,
                    out,
                    width = 7.5,
                    height = 7.5) {
  Cairo::CairoPDF(out, width = width, height = height)
  print(plotobj)
  dev.off()
}


# And now to do the ad-hoc cluster 12 plot
dir <- 'zinzen'
so_markers <- read.table(paste0(dir, '/so_markers.tsv'), sep = '\t')
analyzeCluster(cluster=12, height=2.5, width=7, markerList = so_markers)


  # Now let's do more because the original plots need some tweaks to fit into my PDF

# avalos cluster 7, 10, 8
dir <- 'avalos'
so_markers <- read.table(paste0(dir, '/so_markers.tsv'), sep = '\t')
analyzeCluster(cluster=7, height=2, width=7, markerList = so_markers)
analyzeCluster(cluster=10, height=3, width=7, markerList = so_markers)
analyzeCluster(cluster=8, height=3.5, width=7, markerList = so_markers)


# allen cluster 108, 126
dir <- 'allen'
so_markers <- read.table(paste0(dir, '/so_markers.tsv'), sep = '\t')
analyzeCluster(cluster=108, height=4, width=7, markerList = so_markers)
analyzeCluster(cluster=126, height=3.5, width=7, markerList = so_markers)


#########################################################################
# Make sure my computer doesn't die from all the RAM this script uses, after it's all over
rm(list = ls())
gc()