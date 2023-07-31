########### Round2Anno ###########


#' Get neighbors
#' 
#' Get neighbors for each cell
#'
#' @param train A matrix whose rows are cells and columns are features, each row represents a cell in train dataset
#' @param test A matrix whose rows are cells and columns are features, each row represents a cell to find neighbor from `train`
#' @param k The number of nearest neighbor cells to find for each cell in `test`
#' @param metric The metric to calculate distance, can be `euclidean`, `cosine`, `manhattan`, and `hamming`
#'
#' @return Returns a list containing the nearest neighbor index and the nearest neighbor distance
#' @export
#'
get_neighbors <- function(train,test,k,dist.metric = 'cosine'){
  KNN <- Seurat::FindNeighbors(object = as.matrix(train),query = as.matrix(test),k.param = (k+1),
                               return.neighbor = T,annoy.metric = dist.metric,verbose = FALSE)
  nn.idx <- KNN@nn.idx
  nn.dist = KNN@nn.dist
  rownames(nn.idx) <- rownames(test)
  rownames(nn.dist) <- rownames(test)
  out <- list(nn.idx = nn.idx,
              nn.dist = nn.dist)
  return(out)
}




#' Weighted KNN
#' 
#' Apply the weighted KNN algorithm to make the predictions
#'
#' @param neighbors Neighbors got from `get_neighbors()`
#' @param train.labels The labels of cells in the train dataset
#' @param type Return type, default is NULL, which means only return labels, `prob` means return the predicted probability
#'
#' @return Returns the predicted labels
#' @export
#'
weighted_knn <- function(neighbors,train.labels,type = NULL){
  labels <- as.character(train.labels)
  nn.idx <- neighbors$nn.idx
  nn.dist <- neighbors$nn.dist

  if(is.null(type)){
    scores_with_names <- sapply(1:nrow(nn.idx),function(row_idx){
      one_nn_idx <- nn.idx[row_idx,]
      one_nn_dist <- nn.dist[row_idx,]
      nn_labels <- labels[one_nn_idx]
      if(one_nn_dist[1] < .001){
        one_nn_idx <- one_nn_idx[-1]
        one_nn_dist <- one_nn_dist[-1]
        nn_labels <- nn_labels[-1]
      } else{
        last_idx <- length(one_nn_idx)
        one_nn_idx <- one_nn_idx[-last_idx]
        one_nn_dist <- one_nn_dist[-last_idx]
        nn_labels <- nn_labels[-last_idx]
      }
      one_nn_weight <- (max(one_nn_dist) - one_nn_dist) / (max(one_nn_dist) - min(one_nn_dist))
      labels_impact <- sapply(unique(nn_labels),function(label){
        sum(one_nn_weight[which(nn_labels == label)])
      })
      pred <- names(labels_impact)[which.max(labels_impact)]
      score <- max(labels_impact)/sum(labels_impact)
      names(score) <- pred
      return(score)
    })
    out <- list(pred = names(scores_with_names),
                score = unname(scores_with_names))
  } else{
    score_mtx <- lapply(1:nrow(nn.idx),function(row_idx){
      one_nn_idx <- nn.idx[row_idx,]
      one_nn_dist <- nn.dist[row_idx,]
      nn_labels <- labels[one_nn_idx]
      if(one_nn_dist[1] < .001){
        one_nn_idx <- one_nn_idx[-1]
        one_nn_dist <- one_nn_dist[-1]
        nn_labels <- nn_labels[-1]
      } else{
        last_idx <- length(one_nn_idx)
        one_nn_idx <- one_nn_idx[-last_idx]
        one_nn_dist <- one_nn_dist[-last_idx]
        nn_labels <- nn_labels[-last_idx]
      }
      one_nn_weight <- (max(one_nn_dist) - one_nn_dist) / (max(one_nn_dist) - min(one_nn_dist))
      labels_impact <- sapply(levels(factor(labels)),function(label){
        sum(one_nn_weight[which(nn_labels == label)])
      })
      score <- labels_impact/sum(labels_impact)
    }) %>% do.call(what = rbind)
    rownames(score_mtx) <- rownames(nn.idx)
    pred <- colnames(score_mtx)[apply(score_mtx,1,which.max)]
    score <- apply(score_mtx,1,max)

    out <- list(pred = pred,
                score = score,
                score_mtx = score_mtx)
  }
  return(out)
}




#' Seed cell cleaning
#' 
#' Clean seed cell candidates
#'
#' @param cell_meta A cell metadata 
#' @param query_nmf_embedding Query meta-program matrix
#' @param dist.metric The metric to calculate distance, can be `euclidean`, `cosine`, `manhattan`, and `hamming`
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






#' Predict the label for each query cell using WKNN algorithm
#' 
#' Predict the label for each query cell using WKNN algorithm
#'
#' @param cell_meta A cell meta data
#' @param query_nmf_embedding Query meta-program matrix
#' @param k The number of nearest neighbor cells to find for each query cells
#' @param dist.metric The metric to calculate distance, can be `euclidean`, `cosine`, `manhattan`, and `hamming`
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










