data.dir <- 'data/GSM4132287_Normal_finalaggr'

for (i in seq_along(along.with = data.dir)) {
  run <- data.dir[i]
  if (!dir.exists(paths = run)) {
    stop("Directory provided does not exist")
  }
  barcode.loc <- file.path(run, "barcodes.tsv")
  gene.loc <- file.path(run, "genes.tsv")
  features.loc <- file.path(run, "features.tsv.gz")
  matrix.loc <- file.path(run, "matrix.mtx")
  pre_ver_3 <- file.exists(gene.loc)
  if (!pre_ver_3) {
    addgz <- function(s) {
      return(paste0(s, ".gz"))
    }
    barcode.loc <- addgz(s = barcode.loc)
    matrix.loc <- addgz(s = matrix.loc)
  }
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!pre_ver_3 && !file.exists(features.loc)) {
    stop("Gene name or features file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
}