#' seurat2harmony() Function
#'
#' This function operate harmony process with seurat object (as input data)
#' @param orig.seurat ,dims, harmony.dims, batch.term, nfeatures
#' @keywords seurat2harmony
#' @return
#' @export

seurat2harmony=function(orig.seurat,dims,harmony.dims,batch.term,nfeatures){

  orig.seurat=orig.seurat %>% NormalizeData() %>% FindVariableFeatures(nfeatures=nfeatures) %>%
    ScaleData() %>% RunPCA(npcs=dims) %>%
    FindNeighbors(dims=1:dims) %>% FindClusters(resolution=1) # %>%
  # RunUMAP(dims = 1:dims) %>% RunTSNE(dims = 1:dims,check_duplicate=F)

  harmony.seurat=orig.seurat
  harmony.seurat=RunHarmony(harmony.seurat,batch.term,max.iter.harmony=30)
  harmony.seurat <- harmony.seurat %>%
    RunUMAP(reduction = "harmony", dims = 1:harmony.dims) %>%
    RunTSNE(reduction = "harmony", dims = 1:harmony.dims) %>%
    FindNeighbors(reduction = "harmony", dims = 1:harmony.dims) %>%
    FindClusters(resolution=1) %>% identity()

  return(harmony.seurat)
}
