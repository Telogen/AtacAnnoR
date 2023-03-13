
#' Combine and Discard
#'
#' Combine cell labels predicted by each reference and discard a fraction of cells with uncertain labels
#' @param merged_labels merged_labels got from `get_merged_labels()`
#' @param reference the index of reference to be based on
#' @param fraction default is 0.05
#'
#' @return Return the index of cells to be discarded.
#' @export
#'
Combine_and_Discard <- function(merged_laels,reference = 1,fraction = 0.05,verbose = T){
  
  length_n_labels_df <- table(merged_laels) %>% sort(T) %>% as.data.frame()
  df <- strsplit(as.character(length_n_labels_df$merged_laels),'=') %>% do.call(what = 'rbind') %>% as.data.frame()
  df$num <- length_n_labels_df$Freq
  
  reference1_cts <- unique(df[,reference])
  abandon_idxes <- c()
  abandon_merged_labels <- c()
  for(reference1_ct in reference1_cts){
    # reference1_ct <- '1Monocytes1'
    # print(reference1_ct)
    reference1_ct_df <- df[which(df[,reference] == reference1_ct),]
    reference1_ct_df
    abandon_idx <- rownames(reference1_ct_df[which(reference1_ct_df$num <= sum(reference1_ct_df$num)*fraction),])
    # abandon_idx <- rownames(reference1_ct_df[which(reference1_ct_df$num <= 10),])
    abandon_idxes <- c(abandon_idxes,abandon_idx)
    ct_abandon_merged_labels <- length_n_labels_df[abandon_idx,]$merged_laels %>% as.character()
    abandon_merged_labels <- c(abandon_merged_labels,ct_abandon_merged_labels)
  }
  
  abandon_all_idxes <- which(merged_laels %in% abandon_merged_labels)
  if(verbose){
    message(paste0('Delete ',length(abandon_idxes),' merged labels and ',length(abandon_all_idxes),' cells'))
  }
  return(abandon_all_idxes)
}
