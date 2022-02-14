#' marker_pheatmap() Function
#'
#' This function draws heatmap with cluster-specific marker genes
#' @param geneset,seuratset,scale.max,scale.min,scale.threshold: input parameters
#' @keywords marker_pheatmap
#' @return
#' @export

marker_pheatmap=function(geneset,seuratset,scale.max,scale.min,scale.threshold){

  geneset=geneset[geneset %in% rownames(seuratset)]
  logcounts=GetAssayData(seuratset, slot = "data", assay = "RNA"); logcounts=logcounts[geneset,]

  scaled_mat=t(scale(t(logcounts)))
  scaled_mat=na.omit(scaled_mat)

  max(scaled_mat); min(scaled_mat)
  scaled_mat[scaled_mat > scale.max] <- (scale.max)
  scaled_mat[scaled_mat < (scale.min)]= (scale.min)

  cluster_df=data.frame(row.names=colnames(seuratset), cluster=seuratset$seurat_clusters)

  cluster_num=max(as.numeric(seuratset$seurat_clusters))-1

  dataset=c(); for ( i in 0:cluster_num) {
    dataset=rbind(dataset,rowMeans(scaled_mat[,rownames(cluster_df)[cluster_df$cluster==i]]))
  }; rownames(dataset) <- paste0("Cluster_",c(0:cluster_num))

  dataset=dataset[,(colSums(abs(dataset))>scale.threshold)]

  library(pheatmap)
  mypalette  = colorRampPalette(c("#053061","#2166ac","#f7f7f7", "#d6604d","#9e0142"))(255); range <- max(abs(dataset))
  pheatmap(dataset,color = mypalette,breaks = seq(-range, range, length.out = 255), na_col = "#e0e0e0", show_colnames = T,cluster_cols = T,cluster_rows = T)
}
