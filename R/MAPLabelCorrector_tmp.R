# MAPLabelCorrector_tmp

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



get_multiple_refs_correct_label_0903 <- function(merged_names,cutoff = .5){
  # merged_names <- pred_label_df$merge_name_12345
  merge_name_nmf_pb <- get_pseudo_bulk_mtx(t(query_nmf_embedding),merged_names)
  all_combs_table <- table(merged_names)
  # View(all_combs_table)
  all_combs <- all_combs_table %>% sort(T) %>% names()
  print('all_combs_number')
  print(length(all_combs))

  refnum <- length(strsplit(all_combs[1],'=')[[1]])
  print('refnum')
  print(refnum)

  # 如果有5个标签，将所有尝试限制为最多错3个
  if(refnum <= 4){
    max_edit_distance <- refnum - 1
  } else{
    max_edit_distance <- refnum - 2
  }
  print('max_edit_distance')
  print(max_edit_distance)

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

    # pbmc.ATAC$test <- merged_names
    # plot_highlight_cells(pbmc.ATAC,observe,ident = 'test',label = F) + NoLegend()
    plot_highlight_cells(pbmc.ATAC,'1Naive T cells1=2NK2=3ILC3=4CD8 TCM4=5NK5',ident = 'test',label = F) + NoLegend()
    plot_highlight_cells(pbmc.ATAC,'1Naive T cells1=2CD4n T2=3CD4.Naive3=4CD4 Naive4=5CD4_naive5',ident = 'test',label = F) + NoLegend()
    plot_highlight_cells(pbmc.ATAC,'CD4 Naive',ident = 'true',label = F) + NoLegend()
    plot_highlight_cells(pbmc.ATAC,'gdT',ident = 'true',label = F,pt.size = .2) + NoLegend()

    # MAP
    a <- sapply(to_correct,function(correct_i){
      # correct_i <- to_correct[4]
      # correct_i <- '1Naive T cells1=2CD4n T2=3CD4.Naive3=4CD4 Naive4=5CD4_naive5'
      correct_i

      # P_priori
      P_priori <- all_combs_table[correct_i]/sum(all_combs_table)
      P_priori

      correct_i_idx <- which(merged_names == correct_i)
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
  },mc.cores = 10) %>% do.call(what = c)


  table(out) %>% length() %>% print()
  new_labels <- factor(merged_names,
                       levels = all_combs,
                       labels = out) %>% as.character()
  table(new_labels) %>% sort(T)

  correct_rule <- data.frame(before = all_combs,
                             after = out,
                             changed = all_combs != out)

  new_labels2 <- sapply(new_labels,function(label){
    correct_rule[which(correct_rule$before == label),]$after
  }) %>% unname()
  table(new_labels2) %>% length() %>% print()

  new_labels <- new_labels2
  pbmc.ATAC$new_labels <- new_labels
  pbmc.ATAC$new_label_idx <- factor(new_labels,
                                    levels = table(new_labels) %>% sort(T) %>% names(),
                                    labels = 1:length(table(new_labels)))
  label_idx_df <- data.frame(label = table(new_labels) %>% sort(T) %>% names(),
                             idx = 1:length(table(new_labels)))
  View(label_idx_df)
  # DimPlot(pbmc.ATAC,group.by = 'new_label_idx',label = T) + NoLegend()


  new_1 <- get_new_labels(new_labels,1)
  new_2 <- get_new_labels(new_labels,2)
  new_3 <- get_new_labels(new_labels,3)
  new_4 <- get_new_labels(new_labels,4)
  new_5 <- get_new_labels(new_labels,5)

  out <- list(new_1,new_2,new_3,new_4,new_5)
  return(out)
}



