# data collected from http://scope.aertslab.org/#/f982274f-018a-40a4-94cf-6421f263c4da/Goodwin_Fly_AdultVNC_elife54074.loom/gene
library(tidyverse)
setwd('C:/Users/Edridge/Documents/thesis-edits/analysis/')
# Repetitive code is for cool kids and I'm a cool kid


allenFiles <- list.files(pattern = "Cluster[0-9]{2}.tsv")

allen <- map_df(allenFiles, ~read.table(.x, sep='\t', header=T) %>%
         mutate(cluster = str_remove_all(.x, ".tsv") %>% str_remove_all('Cluster'))
)

dsouza <- read.table('allen_markers.tsv', sep='\t')


allen_soxn <- allen %>% filter(gene=='SoxN') %>% .$cluster %>% unique()
dsouza_soxn <- dsouza %>% filter(gene=='SoxN') %>% .$cluster %>% unique()
allen_d <- allen %>% filter(gene=='D') %>% .$cluster %>% unique()
dsouza_d <- dsouza %>% filter(gene=='D') %>% .$cluster %>% unique()


jaccard <- function(list1, list2) {
  unionSize <- union(list1, list2) %>% unique() %>% length()
  intersectSize <- intersect(list1, list2) %>% unique() %>% length()
  leftSize <- list1 %>% unique() %>% length()
  rightSize <- list2 %>% unique() %>% length()
  leftUnique <- leftSize - intersectSize
  rightUnique <- rightSize - intersectSize
  
  cat(
    paste0(
      '\nLeft total: ', leftSize,
      '\nRight total: ', rightSize,
      '\nTotal elements: ', unionSize,
      '\nShared elements: ', intersectSize,
      '\nJaccard: ',   intersectSize / unionSize,
      '\nLeft unique: ',leftUnique,
      '\nRight unique: ', rightUnique,
      '\nPercent of left list in intersection: ', round(100*intersectSize/leftSize, 4)
    )
  )
}

jaccard(
  allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
)

jaccard(
  allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
)

# For the SoxN-D cluster
jaccard(
  allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster== 108) %>% .$gene %>% unique()
)


###############################################################################
library(readxl)
library(org.Dm.eg.db)
library(rvest)

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



intersectionList <- function(list1, list2) {
  out <- intersect(list1, list2) %>% unique()
  tfs <- unique(out[out %in% tf_list$symbol])
  soxn_bind <- unique(out[out %in% SoxN_bind$Symbol])
  d_bind <- unique(out[out %in% D_bind$Symbol])
  soxn_target <- unique(out[out %in% SoxN_target$Symbol])
  d_target <- unique(out[out %in% D_target$Symbol])  
  
  cat(out, sep='\n')
  cat('_____________________________')
  outputStats <- function(label, subset){
    cat(paste0('\n', label, ': ', length(subset), ' (', round(100*length(subset)/length(out), 4), '%)'))
  }
  cat(paste0('\nTotal size: ', length(out)))
  outputStats('TFs', tfs)
  outputStats('SoxN bound', soxn_bind)
  outputStats('SoxN targets', soxn_target)
  outputStats('D bound', d_bind)
  outputStats('D targets', d_target)
  
  return(out)
}

unionList <- function(list1, list2) {
  out <- union(list1, list2) %>% unique()
  return(out)
}

intersectionList(
  allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
)

intersectionList(
  allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
)

# For the SoxN-D cluster
intersectionList(
  allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster== 108) %>% .$gene %>% unique()
)

########################################################################

lhs <- function(list1, list2) { # The genes unique to the left hand side
  int <- intersect(list1, list2) %>% unique()
  out <- list1[!(list1 %in% int)]
  tfs <- unique(out[out %in% tf_list$symbol])
  soxn_bind <- unique(out[out %in% SoxN_bind$Symbol])
  d_bind <- unique(out[out %in% D_bind$Symbol])
  soxn_target <- unique(out[out %in% SoxN_target$Symbol])
  d_target <- unique(out[out %in% D_target$Symbol]) 
  cat(out, sep='\n')
  cat(paste0('\nNumber of genes: ', length(out)))
  outputStats <- function(label, subset){
    cat(paste0('\n', label, ': ', length(subset), ' (', round(100*length(subset)/length(out), 4), '%)'))
  }
  cat(paste0('\nTotal size: ', length(out)))
  outputStats('TFs', tfs)
  outputStats('SoxN bound', soxn_bind)
  outputStats('SoxN targets', soxn_target)
  outputStats('D bound', d_bind)
  outputStats('D targets', d_target)
  
  return(out)
}

lhs(
  allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
)

lhs(
  allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
)

# For the SoxN-D cluster
lhs(
  allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster== 108) %>% .$gene %>% unique()
)


##############################################################################
# Redo GO analysis for LHS genes, with background being the union of the gene lists
library(topGO)

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

geneSelFn <- function(bgGenes = NA, manualList = NA) {
  genelist <- manualList
  genes <- table(bgGenes) %>% unlist()
  genes[names(genes) %in% genelist] <- T
  genes[!(names(genes) %in% genelist)] <- F
  tblToNumeric(genes, bool = F)
}

selection <- function(bool) {
  return(bool)
}

doAnalysis <-
  function(bgGenes = NA, manualList = NA,  manualListTitle = NA) {
    # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
    bp <-
      annFUN.org(
        whichOnto = "BP",
        feasibleGenes = NULL,
        mapping = "org.Dm.eg.db",
        ID = "symbol"
      )
    
    selectedGenes <- geneSelFn(bgGenes=bgGenes, manualList=manualList)
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
    
    p1title <- manualListTitle
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

bg_soxn <- unionList(
  allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
)
bg_d <- unionList(
  allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
)
bg_both <- unionList(
  allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
  dsouza %>% filter(cluster == 108) %>% .$gene %>% unique()
)



doAnalysis(
  bgGenes = bg_soxn,
  manualList = intersectionList(
    allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
  ),
  manualListTitle = 'SoxN Intersection'
)
doAnalysis(
  bgGenes = bg_soxn,
  manualList = lhs(
    allen %>% filter(cluster %in% allen_soxn) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster %in% dsouza_soxn) %>% .$gene %>% unique()
  ),
  manualListTitle = 'SoxN Allen unique'
)


doAnalysis(
  bgGenes = bg_d,
  manualList = intersectionList(
    allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
  ),
  manualListTitle = 'D Intersection'
)
doAnalysis(
  bgGenes = bg_d,
  manualList = lhs(
    allen %>% filter(cluster %in% allen_d) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster %in% dsouza_d) %>% .$gene %>% unique()
  ),
  manualListTitle = 'D Allen unique'
)


doAnalysis(
  bgGenes = bg_both,
  manualList = intersectionList(
    allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster == 108) %>% .$gene %>% unique()
  ),
  manualListTitle = 'SoxN-D Intersection'
)
doAnalysis(
  bgGenes = bg_both,
  manualList = lhs(
    allen %>% filter(cluster == 25) %>% .$gene %>% unique(),
    dsouza %>% filter(cluster == 108) %>% .$gene %>% unique()
  ),
  manualListTitle = 'SoxN-D Allen unique'
)
