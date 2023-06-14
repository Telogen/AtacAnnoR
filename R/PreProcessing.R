########### PreProcessing ###########
# pre_processing
# mybinarize
# mytfidf
# get_nmf_embedding

#' Pre-process input matrices
#'
#' Pre-process two matrices to remove 0 expressed genes and keep common genes
#'
#' @param ref_mtx a reference gene expression matrix whose rows are genes and columns are cells/samples
#' @param query_mtx query gene activity matrix whose rows are genes and columns are cells
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a list containing two matrices after pre-processing.
#' @export
#'
#'
pre_processing <- function(ref_mtx, query_mtx, verbose = T) {
  ref_mtx <- ref_mtx[which(Matrix::rowSums(ref_mtx) > 0), ]
  query_mtx <- query_mtx[which(Matrix::rowSums(query_mtx) > 0), ]

  used_genes <- intersect(rownames(ref_mtx), rownames(query_mtx))
  if (verbose == T) {
    message(paste0("Keep ", length(used_genes), " intersect genes"))
  }
  ref_mtx <- ref_mtx[used_genes, ]
  query_mtx <- query_mtx[used_genes, ]
  out_list <- list(
    ref_mtx = ref_mtx,
    query_mtx = query_mtx
  )
  return(out_list)
}




mybinarize <- function(object){
  if (inherits(x = object, what = "dgCMatrix")) {
    methods::slot(object = object, name = "x") <-
      rep.int(x = 1, times = length(x = methods::slot(object = object, name = "x")))
  }
  else {
    object[object > 1] <- 1
  }
  return(object)
}


mytfidf <- function (object){
  npeaks <- Matrix::colSums(x = object)
  tf <- Matrix::tcrossprod(x = object, y = Matrix::Diagonal(x = 1/npeaks))
  rsums <- Matrix::rowSums(x = object)
  idf <- ncol(x = object)/rsums
  idf <- log(1 + idf)
  norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
  dimnames(norm.data) <- dimnames(object)
  return(norm.data)
}



#' Get meta-program matrix by NMF
#'
#' Get meta-program matrix by NMF from scATAC-seq peak counts matrix
#'
#' @param peak_counts a query scATAC-seq peak counts whose rows are peaks and columns are cells
#' @param binarize whether to binarize the peak counts matrix, default is TRUE (recommended)
#' @param tfidf whether to do TF-IDF, default is TRUE (recommended)
#' @param normalize whether to do normalization, default is TRUE (recommended)
#' @param nmf_seed the seed set for NMF, default is NULL
#' @param n_factors the number of factors (meta-programs) to get, default is 50
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a dimension reduced matrix whose rows are cells and columns are factors (meta-programs).
#' @export
#'
get_nmf_embedding <- function(peak_counts, binarize = T, tfidf = T, normalize = T,nmf_seed = NULL,
                              n_factors = 50L, verbose = T) {
  if (binarize == T) {
    peak_counts <- mybinarize(peak_counts)
  }

  if (tfidf == T) {
    if (verbose == T) {
      message("Performing TF-IDF normalization...")
    }
    peak_counts <- mytfidf(peak_counts)
  }

  if (normalize == T) {
    peak_counts <- Seurat::NormalizeData(peak_counts, verbose = F)
  }

  if (verbose == T) {
    message("Performing NMF...")
  }
  suppressMessages(library(RcppML))
  ATAC_model <- RcppML::nmf(peak_counts, k = n_factors, nonneg = T, tol = 1e-05, seed = nmf_seed, verbose = F)
  ATAC_emd <- ATAC_model@h
  dimnames(ATAC_emd) <- list(paste0("factor_", 1:n_factors),
                             colnames(peak_counts))
  return(t(ATAC_emd))
}




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
  COR <- pcaPP::cor.fk(as.matrix(pb_counts_selected))
  ord <- stats::hclust(dist(COR, method = "euclidean"), method = "ward.D" )$order
  COR <- COR[ord,ord]
  neighbor_celltypes_list <- apply(COR,1,function(row){
    names(which(sort(row, decreasing = T) > min.cor))
  })
  if (verbose) {
    data <- as.data.frame(as.table(COR))
    data$label <- round((data$Freq*100),0)
    data$label[which(data$label < (min.cor*100))] <- ' '
    p <- ggplot2::ggplot(data, ggplot2::aes(Var1, Var2, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(name = "Similarity\n(×100)",low = "blue", mid = 'yellow', high = "red")  +
      ggplot2::theme_bw() +
      ggplot2::labs(x = " ", y = " ") +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, color = 'black',vjust = 0.5, hjust = 1,size = 10),
                     axis.text.y = ggplot2::element_text(angle = 0,   color = 'black',vjust = 0.5, hjust = 1,size = 10),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.line =  ggplot2::element_blank()) +
      ggplot2::geom_text(ggplot2::aes(label = label),size = 3)
    print(p)
  }
  return(neighbor_celltypes_list)
}


