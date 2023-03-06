########### GetNeighborCTs ###########
# get_pseudo_bulk_mtx
# get_neighbor_celltypes



#' Get the pseudo bulk gene expression matrix
#'
#' Get the pseudo bulk gene expression matrix from a single cell gene expression matrix
#' 
#' @param sc_mtx a single cell gene expression matrix whose rows are features and columns are cells
#' @param labels a vector of cell types whose orders are same as the columns of sc_mtx
#' @param mode the method to aggregate cells, default is \code{'mean'}
#'  \itemize{
#'   \item{'mean': }{get the mean of each cell type}
#'   \item{'sum': }{get the sum of each cell type}
#' }
#'
#' @return Returns a pseudo bulk matrix whose rows are features and columns are cell types.
#' @export
#'
get_pseudo_bulk_mtx <- function(sc_mtx, labels, mode = "mean") {
  Seurat_obj <- SeuratObject::CreateSeuratObject(sc_mtx,meta.data = data.frame(group = labels,row.names = colnames(sc_mtx))) %>% suppressWarnings()
  if (mode == "mean") {
    out_mtx <- Seurat::AverageExpression(Seurat_obj,group.by = 'group',slot = 'counts')[[1]]
  } else{
    out_mtx <- Seurat::AggregateExpression(Seurat_obj,group.by = 'group',slot = 'counts')[[1]]
  }
  return(out_mtx)
}







#' Get neighbor cell types
#'
#' Get neighbor cell types from a pseudo bulk matrix
#'
#' @param sc_count_mtx scRNA-seq counts matrix
#' @param labels scRNA-seq cell labels
#' @param global_markers global markers got from `get_global_markers_sc()`
#' @param min.cor the minimum cutoff of cell type correlation to define 'neighbor cells', default is 0.7
#' @param verbose whether to display messages and plot, default is TRUE
#'
#' @return Returns a list of cell subtypes
#' @export
#'
#' @importFrom pcaPP cor.fk
get_neighbor_celltypes <- function(sc_count_mtx,labels,global_markers,min.cor = 0.7,verbose = T) {
  pb_counts <- get_pseudo_bulk_mtx(sc_count_mtx,labels,mode = 'sum')
  features <- unlist(global_markers,use.names = F) %>% unique()
  pb_counts_selected <- pb_counts[features,]
  COR <- pcaPP::cor.fk(pb_counts_selected)
  neighbor_celltypes_list <- apply(COR,1,function(row){
    names(which(sort(row, decreasing = T) > min.cor))
  })
  if (verbose) {
    ord <- stats::hclust(dist(COR, method = "euclidean"), method = "ward.D" )$order
    data <- as.data.frame(as.table(COR[ord,ord]))
    data$label <- round((data$Freq*100),0)
    data$label[which(data$label < (min.cor*100))] <- ' '
    p <- ggplot2::ggplot(data, ggplot2::aes(Var1, Var2, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(name = "cor",low = "blue", mid = 'yellow', high = "red")  +
      ggplot2::theme_bw() +
      ggplot2::labs(x = " ", y = " ") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line =  ggplot2::element_blank()) +
      ggplot2::geom_text( ggplot2::aes(label = label),size = 3)
    print(p)
  }
  return(neighbor_celltypes_list)
}






