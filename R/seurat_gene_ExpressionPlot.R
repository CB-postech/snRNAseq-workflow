#' seurat_gene_expression_map() Function
#'
#' This function draws gene expression on dimension reduction plot
#' @param seurat.obj,Gene,order,reduction_type,pt.size: input parameters for visualization
#' @keywords seurat_gene_expression_map
#' @return
#' @export

seurat_gene_expression_map <- function(seurat.obj,Gene,order,reduction_type,pt.size) {
  max_value <- max(seurat.obj@assays$RNA@data[Gene, ])
  FeaturePlot(seurat.obj, features = Gene, order = order, reduction = reduction_type, cols = c("#e0e0e0", "#b2182b"), pt.size = pt.size) +
    scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),breaks=c(0, max_value),labels=c(0,ceiling(max_value))) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),plot.title = element_text(Gene),
          legend.position = "none")
}
