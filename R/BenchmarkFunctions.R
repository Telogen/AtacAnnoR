################### BenchmarkFunctions ##################
# get_benchmark
# get_each_recall
# get_each_precision




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



#' Get benchmark indicators
#'
#' @param true_labels a vector of the true labels
#' @param pred_labels a vector of the predicted labels
#'
#' @return Returns a vector of the accuracy, average recall (balanced accuracy), average precision and mean F1 score.
#' @export
#'
get_benchmark <- function(true_labels, pred_labels){
  get_accuracy <- function(true_labels, pred_labels){
    return(length(which(true_labels == pred_labels))/length(true_labels))
  }
  # accuracy
  accuracy <- get_accuracy(true_labels, pred_labels)
  # balanced_accuracy
  df <- data.frame(true = as.character(true_labels), pred = as.character(pred_labels))
  df_li <- split(df,df$true)
  each_group_acc <- sapply(df_li,function(df_i){
    return(get_accuracy(df_i$true,df_i$pred))
  })
  balanced_accuracy <- mean(each_group_acc)
  # average_precision
  df_li2 <- split(df,df$pred)
  each_group_acc2 <- sapply(df_li2,function(df_i){
    return(get_accuracy(df_i$true,df_i$pred))
  })
  average_precision <- mean(each_group_acc2)

  # macro_f1
  class = sort(unique(true_labels))
  tp = NA
  fp = NA
  fn = NA
  for(i in 1:length(class)){
    tp[i] <- sum(pred_labels == class[i] & true_labels == class[i])
    fp[i] <- sum(pred_labels == class[i] & true_labels != class[i])
    fn[i] <- sum(pred_labels != class[i] & true_labels == class[i])
  }
  f1 <- 2*tp/(2*tp+fp+fn)
  names(f1) <- class
  macro_f1 <- mean(f1)

  out <- c(accuracy,balanced_accuracy,average_precision,macro_f1)
  names(out) <- c('accuracy','average_recall','average_precision','macro_f1')
  return(out)
}



#' Get each cell type's recall
#'
#' @param cell_meta the cell metadata
#' @param ident the column name of the predicted labels
#'
#' @return Return a data frame containing each cell type's recall
#' @export
#'
get_each_recall <- function(cell_meta,ident = 'kendall_pred'){
  cell_meta_li <- split(cell_meta,cell_meta$true)
  out <- sapply(cell_meta_li,function(cell_meta_i){
    c(pred_true_cells = length(which(cell_meta_i[,ident] == cell_meta_i$true)),
      all_true_cells = length(cell_meta_i$true),
      recall = length(which(cell_meta_i[,ident] == cell_meta_i$true))/length(cell_meta_i$true))
  }) %>% t() %>% as.data.frame()
  out$pred_true_cells <- as.numeric(out$pred_true_cells)
  out$all_true_cells <- as.numeric(out$all_true_cells)
  out$recall <- as.numeric(out$recall) %>% round(3)
  return(out)
}



#' Get each cell type's precision
#'
#' @param cell_meta the cell metadata
#' @param ident the column name of the predicted labels
#'
#' @return Return a data frame containing each cell type's precision
#' @export
#'
get_each_precision <- function(cell_meta,ident = 'kendall_pred'){
  cell_meta_li <- split(cell_meta,cell_meta[,ident])
  out <- sapply(cell_meta_li,function(cell_meta_i){
    c(true_cells = length(which(cell_meta_i[,ident] == cell_meta_i$true)),
      predicted_cells = length(cell_meta_i$true),
      precision = length(which(cell_meta_i[,ident] == cell_meta_i$true))/length(cell_meta_i$true))
  }) %>% t() %>% as.data.frame()
  out$true_cells <- as.numeric(out$true_cells)
  out$predicted_cells <- as.numeric(out$predicted_cells)
  out$precision <- as.numeric(out$precision) %>% round(3)
  return(out)
}





