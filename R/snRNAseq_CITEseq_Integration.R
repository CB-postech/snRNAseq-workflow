#' snRNAseq_CITEseq_Integration() Function
#'
#' This function adds a and b
#' @param snRNAseq.dir ,CITEseq.dir : data path
#' @keywords snRNAseq_CITEseq_Integration
#' @return
#' @export

snRNAseq_CITEseq_Integration <- function(snRNAseq.dir,CITEseq.dir) {

  library(DropletUtils)
  library(scater)
  library(scRNAseq)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(dplyr)
  library(Matrix)

  set.seed(12345)
  path_snRNAseq=snRNAseq.dir
  snRNAseq.data=Read10X_h5(paste0(path_snRNAseq,"raw_feature_bc_matrix.h5"))
  e.out <- emptyDrops(snRNAseq.data)
  snRNAseq.data.filt=snRNAseq.data[,which(e.out$FDR<0.01)]

  path_HTO=CITEseq.dir
  barcode.file <- paste0(path_HTO, "barcodes.tsv.gz")
  features.file <- paste0(path_HTO, "features.tsv.gz")
  matrix.path <- paste0(path_HTO, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)

  barcode.names = read.delim(barcode.file, header = FALSE, stringsAsFactors = FALSE)
  feature.names = read.delim(features.file, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  colnames(mat) = paste0(colnames(mat),"-1")
  rownames(mat) = feature.names$V1

  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  joint.bcs <- intersect(colnames(snRNAseq.data), colnames(mat))
  snRNAseq.data <- snRNAseq.data[, joint.bcs]
  mat <- as.matrix(mat[, joint.bcs])
  mat.mapped=mat[1:3,] # remove unmapped read row

  # Setup Seurat object
  hashtag <- CreateSeuratObject(counts = snRNAseq.data) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Adding HTO data as an independent assay
  hashtag[["HTO"]] <- CreateAssayObject(counts = mat.mapped)
  hashtag <- NormalizeData(hashtag,assay="HTO",normalization.method="CLR")

  ### Demultiplex cells based on HTO enrichment
  hashtag <- HTODemux(hashtag,assay="HTO")

  # Group cells based on the max HTO signal
  Idents(hashtag) <- "HTO_classification.global"

  # First, we will remove negative cells from the object
  hashtag.singlet <- subset(hashtag,idents="Singlet")

  return (hashtag.singlet)
}

