# Round2Anno
# get_neighbors
# weighted_knn


#
#' get_neighbors
#'
#' @param train todo
#' @param test todo
#' @param k todo
#' @param metric todo
#'
#' @return todo
#' @export
#'
get_neighbors <- function(train,test,k,dist.metric = 'cosine'){
  KNN <- Seurat::FindNeighbors(object = as.matrix(train),query = as.matrix(test),k.param = (k+1),
                               return.neighbor = T,annoy.metric = dist.metric,verbose = F)
  nn.idx <- KNN@nn.idx
  nn.dist = KNN@nn.dist
  rownames(nn.idx) <- rownames(test)
  rownames(nn.dist) <- rownames(test)
  out <- list(nn.idx = nn.idx,
              nn.dist = nn.dist)
  return(out)
}




#' weighted_knn
#'
#' @param neighbors todo
#' @param train.labels todo
#' @param type NULL or 'prob'
#'
#' @return todo
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



#' get_final_seed_barcodes
#'
#' @param cell.meta todo
#' @param candidate.seed.barcodes todo
#' @param query.nmf.embedding todo
#'
#' @return todo
#' @export
#'
#' @examples todo
get_final_seed_barcodes <- function(cell.meta,candidate.seed.barcodes,query.nmf.embedding, k = 10){
  seed_meta <- cell.meta[candidate.seed.barcodes,]
  seed_neighbors <- get_neighbors(query.nmf.embedding[candidate.seed.barcodes,],
                                  query.nmf.embedding[candidate.seed.barcodes,],k = 10)
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
  return(final_seed_barcodes)
}










