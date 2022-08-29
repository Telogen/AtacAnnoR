########### PreProcessing ###########
# pre_processing
# get_nmf_embedding
# get_pseudo_bulk_mtx
# get_pb_markers




#' Pre-process input matrices
#'
#' Pre-process two matrices to remove 0 expressed genes and keep common genes
#'
#' @param ref.mtx a reference gene expression matrix whose rows are genes and columns are cells/samples
#' @param query.mtx query gene activity matrix whose rows are genes and columns are cells
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a list containing two matrices after pre-processing.
#' @export
#'
#'
pre_processing <- function(ref.mtx, query.mtx, verbose = T) {
  ref.mtx <- ref.mtx[which(Matrix::rowSums(ref.mtx) > 0), ]
  query.mtx <- query.mtx[which(Matrix::rowSums(query.mtx) > 0), ]

  used_genes <- intersect(rownames(ref.mtx), rownames(query.mtx))
  if (verbose == T) {
    message(paste0("Keep ", length(used_genes), " intersect genes"))
  }
  ref.mtx <- ref.mtx[used_genes, ]
  query.mtx <- query.mtx[used_genes, ]
  out_list <- list(
    ref.mtx = ref.mtx,
    query.mtx = query.mtx
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
  norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
  dimnames(norm.data) <- dimnames(object)
  return(norm.data)
}



#' Get dimension reduced cell embedding by NMF
#'
#' Get dimension reduced cell embedding by NMF from scATAC-seq peak counts
#'
#' @param peak.counts a query scATAC-seq peak counts whose rows are peaks and columns are cells
#' @param binarize whether to binarize the peak counts matrix, default is TRUE (recommended)
#' @param tfidf whether to do TF-IDF, default is TRUE (recommended)
#' @param normalize whether to do normalization, default is TRUE (recommended)
#' @param nmf.seed the seed set for \code{RcppML::nmf}, default is NULL
#' @param n.factors the number of factors (meta-programs) to get, default is 50
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a dimension reduced matrix whose rows are cells and columns are
#' factors (meta-programs).
#' @export
#'
get_nmf_embedding <- function(peak.counts, binarize = T, tfidf = T, normalize = T,nmf.seed = NULL,
                              n.factors = 50L, verbose = T) {
  if (binarize == T) {
    peak.counts <- mybinarize(peak.counts)
  }

  if (tfidf == T) {
    if (verbose == T) {
      message("Performing TF-IDF normalization...")
    }
    peak.counts <- mytfidf(peak.counts)
  }

  if (normalize == T) {
    peak.counts <- Seurat::NormalizeData(peak.counts, verbose = F)
  }

  if (verbose == T) {
    message("Performing NMF...")
  }
  suppressMessages(library(RcppML))
  ATAC_model <- RcppML::nmf(peak.counts, k = n.factors, nonneg = T, tol = 1e-05, seed = nmf.seed, verbose = F)
  ATAC_emd <- ATAC_model$h
  dimnames(ATAC_emd) <- list(paste0("factor_", 1:n.factors),
                             colnames(peak.counts))
  return(t(ATAC_emd))
}




#' Get the pseudo bulk matrix
#'
#' Get the pseudo bulk matrix from a single cell matrix
#'
#' @param sc.mtx a single cell matrix whose rows are features and columns are cells
#' @param cell.type a vector of cell types whose orders are same as the columns of sc.mtx
#' @param mode the method to aggragate cells, default is \code{'mean'}
#'  \itemize{
#'   \item{'mean': }{get the mean of each cell type}
#'   \item{'sum': }{get the sum of each cell type}
#'   \item{'median': }{get the median of each cell type}
#' }
#'
#' @return Returns a pseudo bulk matrix whose rows are features and columns are cell types.
#' @export
#'
get_pseudo_bulk_mtx <- function(sc.mtx, cell.type, mode = "mean") {
  cell.type <- as.character(cell.type)
  sapply(unique(cell.type), function(ct) {
    ct_mtx <- sc.mtx[, which(cell.type == ct)]
    if (class(ct_mtx)[1] != "numeric") {
      if (mode == "mean") {
        Matrix::rowSums(ct_mtx) / ncol(ct_mtx) %>% return()
      } else if (mode == "sum") {
        Matrix::rowSums(ct_mtx) %>% return()
      } else if (mode == "median") {
        apply(ct_mtx, 1, stats::median) %>% return()
      }
    } else {
      ct_mtx %>% return()
    }
  }) %>%
    return()
}



#' Get markers from a pseudo bulk matrix
#'
#' Get markers from a pseudo bulk matrix
#'
#' @param pb.mtx a pseudo bulk matrix whose rows are features and columns are cell types
#' @param ct.to.find.mk a vector of cell types to find markers, default is NULL (all cell types)
#' @param HEG.prop the proportion of highly expressed genes to first find, default is 0.2,
#' if the matrix is a real bulk matrix, this parameter is recommended to set as 0.8
#' @param each.ct.top the number of top marker genes to find for each cell type
#' @param return.type the type to return, default is \code{'df'}
#'  \itemize{
#'   \item{df: }{a dataframe whose rows are marker genes and columns are cell types}
#'   \item{vector: }{a vector contains marker genes of all cell types}
#' }
#'
#' @return If \code{return.type = 'df'}, returns a data frame whose rows are marker genes and columns are cell types.
#'  If \code{return.type = 'vector'}, returns a vector containing marker genes of all cell types.
#' @export
#'
get_pb_markers <- function(pb.mtx, ct.to.find.mk = NULL, HEG.prop = 0.2, each.ct.top = 200,min.cor = NULL,return.cor.mtx = FALSE) {
  feature.use <- apply(pb.mtx, 2, function(col) {
    sort(col, T)[1:ceiling(nrow(pb.mtx) * HEG.prop)] %>%
      names()
  }) %>% as.vector() %>% unique()
  pb.mtx.use <- pb.mtx[feature.use, ]
  ideal_gene_mtx <- diag(ncol(pb.mtx.use))
  colnames(ideal_gene_mtx) <- colnames(pb.mtx.use)
  if (!is.null(ct.to.find.mk)) {
    ideal_gene_mtx <- ideal_gene_mtx[, ct.to.find.mk]
  }
  cos_mtx <- RcppML::cosine(t(pb.mtx.use), ideal_gene_mtx)
  rownames(cos_mtx) <- rownames(pb.mtx.use)
  colnames(cos_mtx) <- colnames(ideal_gene_mtx)
  if(is.null(min.cor)){
    out <- apply(cos_mtx, 2, function(col) {
      sort(col, T)[1:each.ct.top] %>% names()
    }) %>% as.data.frame()
  } else{
    out <- apply(cos_mtx,2,function(col){
      rownames(cos_mtx)[which(col > min.cor)]
    })
  }
  if(return.cor.mtx){
    return(list(markers = out,
                cor_mtx = as.data.frame(cos_mtx)))
  } else{
    return(out)
  }
}
