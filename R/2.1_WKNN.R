########### WKNN ###########
# get_neighbors
# weighted_knn


#' get_neighbors
#'
#' @param train todo
#' @param test todo
#' @param k todo
#' @param metric euclidean, cosine, manhattan, and hamming
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











