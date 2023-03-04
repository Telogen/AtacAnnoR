########### WKNNPredict ###########
# WKNN_predict

WKNN_predict <- function(cell_meta,query_nmf_embedding,k = 10,dist.metric = 'cosine'){
  final_seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  final_seed_barcodes <- rownames(final_seed_meta)
  
  neighbors <- get_neighbors(query_nmf_embedding[final_seed_barcodes,],query_nmf_embedding,k = k,dist.metric = dist.metric)
  final_pred <- weighted_knn(neighbors,final_seed_meta$kendall_pred)
  cell_meta$final_pred <- final_pred$pred
  return(cell_meta)
}


