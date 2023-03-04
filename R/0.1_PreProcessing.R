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
#' @examples
#' pre_processing_mtxs <- pre_processing(ref_mtx = REF_MTX, query_mtx = GA_MTX)
#' ref_mtx <- pre_processing_mtxs$ref_mtx
#' query_mtx <- pre_processing_mtxs$query_mtx
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



#' Get dimension reduced cell embedding by NMF
#'
#' Get dimension reduced cell embedding by NMF from scATAC-seq peak counts
#'
#' @param peak_counts a query scATAC-seq peak counts whose rows are peaks and columns are cells
#' @param binarize whether to binarize the peak counts matrix, default is TRUE (recommended)
#' @param tfidf whether to do TF-IDF, default is TRUE (recommended)
#' @param normalize whether to do normalization, default is TRUE (recommended)
#' @param nmf_seed the seed set for \code{RcppML::nmf}, default is NULL
#' @param n_factors the number of factors (meta-programs) to get, default is 50
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a dimension reduced matrix whose rows are cells and columns are
#' factors (meta-programs).
#' @export
#'
#' @examples
#' query_nmf_embedding <- get_nmf_embedding(PEAK_COUNTS, n_factors = 50)
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
  ATAC_emd <- ATAC_model$h
  dimnames(ATAC_emd) <- list(paste0("factor_", 1:n_factors),
                             colnames(peak_counts))
  return(t(ATAC_emd))
}




