# MAPLabelCorrector

get_edit_dist_n_combs <- function(one_comb,all_combs,dist_n = 1){
  labels <- strsplit(one_comb,'=')[[1]]
  regular_tmp <- combn(labels,(length(labels) - dist_n),simplify = F)
  out <- sapply(regular_tmp,function(one_regular_tmp){
    one_regular_exp <- paste(one_regular_tmp,collapse = '.*') %>%
      gsub(pattern = '\\+',replacement = '\\.\\*')
    all_combs[grep(one_regular_exp,all_combs)]
  }) %>% unlist() %>% unique()
  return(out)
}




#' get merged labels
#'
#' @param pred_df 
#'
#' @return Return the merged label
#' @export
#'
get_merged_labels <- function(pred_df){
  pred_df_with_num <- sapply(1:ncol(pred_df),function(i){
    paste0(i,pred_df[[i]],i)
  }) %>% as.data.frame()
  
  merged_labels <- apply(pred_df_with_num,1,function(row){
    paste(row,collapse = '=')
  })
  return(merged_labels)
}


#' correct labels
#'
#' @param merged_labels merged_labels got from `get_merged_labels()`
#' @param query_nmf_embedding the nmf_embedding of query data
#' @param cutoff 0.6 by default, smaller value means more strict correction
#' @param max_edit_distance max edit distance, by default is the number of references-1 
#' @param threads number of threads
#' @param verbose verbose
#'
#' @return Return a list containing the information of the corrected labels  
#' @export
#'
correct_labels <- function(merged_labels,query_nmf_embedding,cutoff = .6,max_edit_distance = NULL,threads = 10, verbose = T){
  # merged_labels <- pred_label_df$merge_name_12345
  merge_name_nmf_pb <- get_pseudo_bulk_mtx(t(query_nmf_embedding),merged_labels)
  all_combs_table <- table(merged_labels)
  # View(all_combs_table)
  all_combs <- all_combs_table %>% sort(T) %>% names()
  refnum <- length(strsplit(all_combs[1],'=')[[1]])
  
  if(is.null(max_edit_distance)){
    # max_edit_distance <- ceiling(refnum/2)
    max_edit_distance <- refnum - 1
  }
  
  if(verbose){
    message(paste0('Cluster number before correction: ',length(all_combs)))
    message(paste0('Reference number: ',refnum))
    message(paste0('Max edit distance: ',max_edit_distance))
  }
  
  out <- pbmcapply::pbmclapply(all_combs,function(one_comb){
    # one_comb <- all_combs[151]
    # one_comb <- '1Naive T cells1=2CD4n T2=3CD4.Th3=4Treg4=5TFH5'
    one_comb
    to_correct <- get_edit_dist_n_combs(one_comb,all_combs,max_edit_distance)
    to_correct
    all_combs_table[to_correct] %>% sort(T) %>% head()
    all_combs_table[one_comb]
    observe <- one_comb
    observe

    # SeuratObj_ATAC$test <- merged_labels
    # plot_highlight_cells(SeuratObj_ATAC,observe,ident = 'test',label = F) + NoLegend()
    # plot_highlight_cells(SeuratObj_ATAC,'1Naive T cells1=2NK2=3ILC3=4CD8 TCM4=5NK5',ident = 'test',label = F) + NoLegend()
    # plot_highlight_cells(SeuratObj_ATAC,'1Naive T cells1=2CD4n T2=3CD4.Naive3=4CD4 Naive4=5CD4_naive5',ident = 'test',label = F) + NoLegend()
    # plot_highlight_cells(SeuratObj_ATAC,'CD4 Naive',ident = 'true',label = F) + NoLegend()
    # plot_highlight_cells(SeuratObj_ATAC,'gdT',ident = 'true',label = F,pt.size = .2) + NoLegend()

    # MAP
    a <- sapply(to_correct,function(correct_i){
      # correct_i <- to_correct[4]
      # correct_i <- '1Naive T cells1=2CD4n T2=3CD4.Naive3=4CD4 Naive4=5CD4_naive5'
      correct_i

      # P_priori
      P_priori <- all_combs_table[correct_i]/sum(all_combs_table)
      P_priori

      correct_i_idx <- which(merged_labels == correct_i)
      correct_i_idx
      correct_i_nmf_embedding <- query_nmf_embedding[correct_i_idx,]
      observe_pb_embedding <- merge_name_nmf_pb[,observe]
      if(length(correct_i_idx) == 1){
        Coss <- RcppML::cosine(correct_i_nmf_embedding,observe_pb_embedding)
      } else{
        Coss <- RcppML::cosine(t(correct_i_nmf_embedding),observe_pb_embedding) # %>% mean()
      }
      # P_likelihood <- RcppML::cosine(merge_name_nmf_pb[,observe],merge_name_nmf_pb[,correct_i])
      P_likelihood <- length(which(Coss > cutoff))/length(Coss)
      P_likelihood
      out <- P_priori*P_likelihood

      return(out)
      # out
      # P_likelihood
      # Coss
      # P_priori
    })
    names(a) <- to_correct
    sort(a,T) %>% head()
    this_out <- sort(a,T)[1] %>% names()
    return(this_out)
  },mc.cores = threads) %>% do.call(what = c)
  
  correct_rule2 <- data.frame(before = all_combs,
                              freq = all_combs_table %>% sort(T) %>% as.numeric(),
                              after = out) %>% 
    mutate(changed = before != after)
  
  # make correct_rule
  correct_rule <- out
  names(correct_rule) <- all_combs
  
  # correct labels iteratively using correct_rule: get new_combs
  before_num <- length(unique(all_combs))
  diff <- before_num
  new_combs <- all_combs
  while(diff > 0){
    new_combs <- correct_rule[new_combs]
    after_num <- length(unique(new_combs))
    # print(after_num)
    diff <- before_num - after_num
    # print(diff)
    before_num <- after_num
    # print('=====')
  }
  print(paste0('Cluster number after correction: ',after_num))
  
  # get new labels
  new_labels <- factor(merged_labels,
                       levels = all_combs,
                       labels = new_combs) %>% as.character()
  
  # corrected label indexes
  new_label_idx <- factor(new_labels,
                          levels = table(new_labels) %>% sort(T) %>% names(),
                          labels = 1:length(table(new_labels)))
  # corrected label and indexes df
  label_idx_df <- data.frame(label = table(new_labels) %>% sort(T) %>% names(),
                             idx = 1:length(table(new_labels)),
                             num = (table(new_labels) %>% sort(T) %>% as.numeric()))
  label_idx_df
  # each corrected labels
  get_new_labels <- function(merged_labels,i){
    merged_labels_list <- strsplit(merged_labels,'=')
    sapply(merged_labels_list,function(x){substr(x[i],2,nchar(x[i])-1)})
  }
  
  corrected_df <- data.frame(row.names = 1:length(new_labels))
  for(i in 1:refnum){
    corrected_df[[i]] <- get_new_labels(new_labels,i)
  }

  out <- list(corrected_merged_labels = new_labels,
              corrected_merged_labels_idx = new_label_idx,
              label_idx_df = label_idx_df,
              each_corrected_df = corrected_df,
              correct_rule = correct_rule2)
  return(out)
}



