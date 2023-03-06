########### WKNNPredict ###########
# WKNN_predict

#' Predict the label for each query cell using WKNN algorithm
#' 
#' Predict the label for each query cell using WKNN algorithm
#'
#' @param cell_meta a cell meta data
#' @param query_nmf_embedding query meta-program matrix
#' @param k the number of nearest neighbor cells to find for each query cells
#' @param dist.metric the metric to calculate distance, can be `euclidean`, `cosine`, `manhattan`, and `hamming`, default is `cosine`
#'
#' @return Returns a new cell meta data with a column `final_pred`
#' @export
#'
WKNN_predict <- function(cell_meta,query_nmf_embedding,k = 10,dist.metric = 'cosine'){
  final_seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  final_seed_barcodes <- rownames(final_seed_meta)
  
  neighbors <- get_neighbors(query_nmf_embedding[final_seed_barcodes,],query_nmf_embedding,k = k,dist.metric = dist.metric)
  final_pred <- weighted_knn(neighbors,final_seed_meta$kendall_pred)
  cell_meta$final_pred <- final_pred$pred
  return(cell_meta)
}


