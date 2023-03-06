########### SeedCleaning ###########
# seed_cleaning


#' Seed cell cleaning
#' 
#' Clean seed cell candidates
#'
#' @param cell_meta a cell metadata 
#' @param query_nmf_embedding query meta-program matrix
#' @param dist.metric the metric to calculate distance, can be `euclidean`, `cosine`, `manhattan`, and `hamming`, default is `cosine`
#'
#' @return Returns a new cell meta data with a column `is_seed`
#' @export
#'
seed_cleaning <- function(cell_meta,query_nmf_embedding,dist.metric = 'cosine'){
  
  seed_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
  seed_barcodes <- rownames(seed_meta)
  
  seed_neighbors <- get_neighbors(query_nmf_embedding[seed_barcodes,],
                                  query_nmf_embedding[seed_barcodes,],k = 10,dist.metric = dist.metric)
  pred2 <- weighted_knn(seed_neighbors,seed_meta$kendall_pred)
  seed_meta$pred_score <- pred2$score
  final_seed_meta_tmp <- seed_meta[which((seed_meta$kendall_pred == pred2$pred)),]
  final_seed_meta_li <- split(final_seed_meta_tmp,final_seed_meta_tmp$kendall_pred)
  final_seed_barcodes <- lapply(final_seed_meta_li,function(final_seed_meta_i){
    if(nrow(final_seed_meta_i) > 200){
      keep_barcode <- rownames(final_seed_meta_i)[which(final_seed_meta_i$pred_score > .9)]
    } else{
      keep_barcode <- rownames(final_seed_meta_i)
    }
    return(keep_barcode)
  }) %>% unlist(use.names = F)
  
  cell_meta$is_seed <- FALSE
  cell_meta[final_seed_barcodes,]$is_seed <- TRUE
  
  return(cell_meta)
}



