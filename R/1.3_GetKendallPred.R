########### GetKendallPred ###########
# get_cor_mtx
# get_kendall_pred
# get_cell_meta_with_true





#' Get correlation matrix
#'
#' Get correlation matrix of each query cell to each reference cell type using Kendall correlation coefficient
#'
#' @param sc_count_mtx reference gene counts matrix
#' @param labels reference cell type labels
#' @param query_mtx query_mtx query gene activity matrix
#' @param global_markers global markers got from `get_global_markers_sc()` or `get_global_markers_bulk()`
#' @param query_nmf_embedding query meta-program matrix
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a correlation matrix whose rows are query cells and columns are reference cell types
#' @export
#'
get_cor_mtx <- function(sc_count_mtx,labels,query_mtx,global_markers,query_nmf_embedding,threads = 10,verbose = T) {
  # ref
  if (verbose) {
    message("Processing reference...")
  }
  ref_selected_features <- unlist(global_markers) %>% unique()
  pb_ref_counts <- get_pseudo_bulk_mtx(sc_count_mtx,labels,mode = 'sum')
  # query
  if (verbose) {
    message("Processing query...")
  }
  query_clusters <- (Seurat::FindNeighbors(query_nmf_embedding, verbose = F)$snn %>%
                       Seurat::FindClusters(verbose = F))[,1] %>% as.character()
  query_global_markers <- get_global_markers_sc(sc_counts_mtx = query_mtx,labels = query_clusters)
  query_selected_features <- unlist(query_global_markers,use.names = F) %>% unique()
  # query_selected_features <- rownames(query_mtx)
  
  # intersect_genes
  intersect_genes <- intersect(ref_selected_features, query_selected_features)
  feature_selected_pb_ref <- pb_ref_counts[intersect_genes, ]
  feature_selected_query_mtx <- query_mtx[intersect_genes, ]
  
  if (verbose) {
    message(paste0("Intersect genes number: ", length(intersect_genes)))
    message(paste0("Computing Kendall correlation coefficient using ", threads, " threads..."))
  }
  # correlation
  celltypes <- colnames(feature_selected_pb_ref)
  out <- pbmcapply::pbmclapply(celltypes, function(celltype) {
    celltype_exp <- feature_selected_pb_ref[, celltype]
    celltype_cor <- apply(feature_selected_query_mtx, 2, pcaPP::cor.fk, y = celltype_exp)
    return(celltype_cor)
  }, mc.cores = threads) %>% do.call(what = 'cbind')
  colnames(out) <- celltypes
  out <- as.data.frame(out)
  return(out)
}





#' Get first round annotation labels
#'
#' Get first round annotation labels from the correlation matrix
#'
#' @param cor_mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a cell metadata whose rows are query cells and columns are first round labels.
#' @export
#'
#'
get_kendall_pred <- function(cor_mtx) {
  cell_meta <- data.frame(row.names = rownames(cor_mtx))
  cell_meta$kendall_pred <- colnames(cor_mtx)[apply(cor_mtx, 1, function(row) {
    ifelse(any(is.na(row)), c(1), which.max(row))
  })]
  return(cell_meta)
}


#' Add the true labels to cell metadata
#' 
#' Add the true labels to cell metadata
#'
#' @param cell_meta a cell metadata
#' @param true_labels a vector of true cell labels
#' @param cor_mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a new cell metadata with true labels.
#' @export
#'
#'
get_cell_meta_with_true <- function(cell_meta,true_labels,cor_mtx = NULL) {
  if (!is.null(cor_mtx)) {
    get_rank <- function(row, names, ct.num, whose.rank) {
      names(row) <- names
      row_cor <- as.numeric(row[1:ct.num])
      names(row_cor) <- names[1:ct.num]
      sorted_row <- sort(row_cor, decreasing = T)
      rank <- which(names(sorted_row) == row[whose.rank])
      return(rank)
    }
    ct_num <- ncol(cor_mtx)
    cor_mtx$true <- true_labels
    cell_meta$true_rank <- apply(cor_mtx, 1, get_rank, names = colnames(cor_mtx), ct_num, whose.rank = "true")
  }
  
  cell_meta$true <- true_labels
  if ('kendall_pred' %in% colnames(cell_meta)){
    cell_meta$kendall_pred_booltrue <- (cell_meta$kendall_pred == true_labels)
  }
  return(cell_meta)
}





