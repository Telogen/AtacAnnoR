################### UsefulFunctions ##################
# get_benchmark
# get_each_accuracy
# plot_confusion_matrix
# plot_highlight_cells
# plot_cell_distribution
# plot_seeds
# plot_pred_scores
# get_signature
# get_MAESTRO_signatures


#' Get benchmark indicators
#'
#' @param true.labels a vector of the true labels
#' @param pred.labels a vector of the predicted labels
#'
#' @return Returns a vector of the accuracy, average recall (balanced accuracy), average precision and mean F1 score.
#' @export
#'
get_benchmark <- function(true.labels, pred.labels){
  get_accuracy <- function(true.labels, pred.labels){
    return(length(which(true.labels == pred.labels))/length(true.labels))
  }
  # accuracy
  accuracy <- get_accuracy(true.labels, pred.labels)
  # balanced_accuracy
  df <- data.frame(true = as.character(true.labels), pred = as.character(pred.labels))
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
  class = sort(unique(true.labels))
  tp = NA
  fp = NA
  fn = NA
  for(i in 1:length(class)){
    tp[i] <- sum(pred.labels == class[i] & true.labels == class[i])
    fp[i] <- sum(pred.labels == class[i] & true.labels != class[i])
    fn[i] <- sum(pred.labels != class[i] & true.labels == class[i])
  }
  f1 <- 2*tp/(2*tp+fp+fn)
  names(f1) <- class
  macro_f1 <- mean(f1)

  out <- c(accuracy,balanced_accuracy,average_precision,macro_f1)
  names(out) <- c('accuracy','average_recall','average_precision','macro_f1')
  return(out)
}



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




plot_confusion_matrix <- function(cell_meta, which.group, mode = 1,title = NULL){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(cowplot) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()

  bench <- get_benchmark(cell_meta$true,cell_meta[,which(colnames(cell_meta) == which.group)])
  bench <- round(bench*100,2)
  if(is.null(title)){
    title <- mode
  }

  if(title == 1){
    myggtitle <- ggtitle(paste0(which.group,' average recall: ',bench[2],'%'))
  } else if(title == 2){
    myggtitle <- ggtitle(paste0(which.group,' average precision: ',bench[3],'%'))
  } else if(title == 3){
    myggtitle <- ggtitle(paste0(which.group,' accuracy: ',bench[1],'%'))
  }

  # calculate recall
  if (mode == 1){
    predictions <- table(cell_meta$true, cell_meta[,which(colnames(cell_meta) == which.group)])
    predictions <- predictions/rowSums(predictions)
    predictions <- as.data.frame(predictions)
    common_levels <- intersect(levels(predictions$Var1),levels(predictions$Var2))
    predictions$Var1 <- factor(predictions$Var1,
                               levels = c(common_levels,setdiff(levels(predictions$Var1),common_levels)))
    predictions$Var2 <- factor(predictions$Var2,
                               levels = c(common_levels,setdiff(levels(predictions$Var2),common_levels)))
    predictions$freq_morethan5 <- round(predictions$Freq,2)
    predictions$freq_morethan5[which(predictions$freq_morethan5 < .05)] <- ' '

    g <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
      geom_tile() +
      myggtitle +
      scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025",limits = c(0,1))  +
      xlab("True cell type") +
      ylab("Predicted cell type") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_text(aes(label = freq_morethan5),size = 3)

  }

  # calculate precision
  if (mode == 2){
    predictions <- table(cell_meta$true, cell_meta[,which(colnames(cell_meta) == which.group)])
    predictions <- t(predictions)/colSums(predictions)
    predictions <- as.data.frame(predictions)
    common_levels <- intersect(levels(predictions$Var1),levels(predictions$Var2))
    predictions$Var1 <- factor(predictions$Var1,
                               levels = c(common_levels,setdiff(levels(predictions$Var1),common_levels)))
    predictions$Var2 <- factor(predictions$Var2,
                               levels = c(common_levels,setdiff(levels(predictions$Var2),common_levels)))
    predictions$freq_morethan5 <- round(predictions$Freq,2)
    predictions$freq_morethan5[which(predictions$freq_morethan5 < .05)] <- ' '

    g <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
      geom_tile() +
      myggtitle +
      scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025",limits = c(0,1))  +
      xlab("Predicted cell type") +
      ylab("True cell type") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_text(aes(label = freq_morethan5),size = 3)

  }

  # number
  if (mode == 3){
    predictions <- table(cell_meta$true, cell_meta[,which(colnames(cell_meta) == which.group)])
    predictions <- as.data.frame(predictions)
    common_levels <- intersect(levels(predictions$Var1),levels(predictions$Var2))
    predictions$Var1 <- factor(predictions$Var1,
                               levels = c(common_levels,setdiff(levels(predictions$Var1),common_levels)))
    predictions$Var2 <- factor(predictions$Var2,
                               levels = c(common_levels,setdiff(levels(predictions$Var2),common_levels)))
    g <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
      geom_tile() +
      myggtitle +
      scale_fill_gradient(name = "Number of cells",low = "#ffffc8", high = "#7d0025")  +
      ylab("Predicted cell type") +
      xlab("True cell type") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_text(aes(label = Freq),size = 3)

  }

  return(g)
}



plot_highlight_cells <- function(Seurat.object, celltype, ident, label = T, pt.size = 1,reduction = DefaultDimReduc(Seurat.object)){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()

  Idents(Seurat.object) <- ident
  cells.highlight <- list()
  for (i in celltype){
    cells.highlight[[`i`]] <- WhichCells(Seurat.object, idents = i)
  }
  p <- DimPlot(Seurat.object, label = label, group.by = ident, cells.highlight = cells.highlight, reduction = reduction,
               cols.highlight = c("red", "orange",'yellow','green','skyblue','pink2','purple','darkgreen','blue'),
               sizes.highlight = pt.size,repel = T) + ggtitle(ident)
  return(p)
}




plot_cell_distribution <- function(cell_meta,celltypes.to.plot = NULL,group.by = 'kendall_pred',
                                   x = 'NMSS',y = 'GMSS',shape = 'is_seed'){
  library(aplot) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  if(is.null(cell_meta$is_seed)){
    cell_meta$is_seed <- FALSE
  }
  name_booltrue <- paste0(group.by,'_booltrue')
  cell_meta[[name_booltrue]] <- (cell_meta[,group.by] == cell_meta[,'true'])

  cell_meta_li <- split(cell_meta,cell_meta[,group.by])

  if (is.null(celltypes.to.plot)){
    cell_meta_li <- cell_meta_li
  } else{
    cell_meta_li <- cell_meta_li[celltypes.to.plot]
  }

  for(i in cell_meta_li){
    DATA <- i

    p1 <- ggplot(data = DATA) +
      geom_point(aes_string(x = x,y = y, color = name_booltrue,shape = shape)) +
      scale_size(range = c(0.1, 3))
    # + geom_smooth(aes(x = NMSS,y = NMFC),method = 'lm',formula = 'y ~ x',size = .5)


    p2 <- ggplot(DATA,aes_string(x = x,color = name_booltrue)) +
      geom_density(aes(y = after_stat(count),fill = .data[[name_booltrue]],alpha = 0.2)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void()+
      theme(axis.title = element_blank(),
            axis.text = element_blank()) +
      NoLegend() +
      ggtitle(paste0(DATA[1,group.by],
                     ' T:',table(DATA[,name_booltrue])['TRUE'],
                     ' F:',table(DATA[,name_booltrue])['FALSE'])) +
      theme(plot.title = element_text(hjust = 0)) +
      theme(plot.title = element_text(vjust = -3))

    p3 <- ggplot(DATA,aes_string(x = y,color = name_booltrue)) +
      geom_density(aes(y = after_stat(count),fill = .data[[name_booltrue]],alpha = 0.2)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      coord_flip() + NoLegend()

    p <- p1 %>% insert_top(p2,height = 0.2) %>% insert_right(p3,0.2)

    print(p)
  }
}




plot_seeds <- function(Seurat.object,cell_meta,ident,pt.size = 0.5,label = F,reduction = DefaultDimReduc(Seurat.object)){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()

  seed.cell_meta <- cell_meta[which(cell_meta$is_seed == TRUE),]
  seed.cell_meta$booltrue <- seed.cell_meta[[ident]] == seed.cell_meta$true
  cells <- list(true_seeds = rownames(seed.cell_meta[which(seed.cell_meta$booltrue == 'TRUE'),]),
                false_seeds = rownames(seed.cell_meta[which(seed.cell_meta$booltrue == 'FALSE'),]))
  p <- DimPlot(Seurat.object,
               cells.highlight = cells,
               group.by = "true",
               reduction = reduction,
               sizes.highlight = pt.size,
               label = label,
               cols.highlight = c('red','blue')) +
    ggtitle(paste0('Seed cell accuracy: ',
                   round(get_benchmark(seed.cell_meta$true,seed.cell_meta[[ident]])[1],3)))
  return(p)
}




plot_pred_scores <- function(cell_meta,mode = 1){
  cell_meta$annotation_correct <- factor(cell_meta$true == cell_meta$final_pred,levels = c(TRUE,FALSE))
  if (mode == 1){
    p1 <- ggplot(cell_meta, aes(pred_score, fill = annotation_correct, colour = annotation_correct)) +
      geom_density(position = "fill",aes(y = after_stat(count),alpha = 0.5)) +
      scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      xlab("Prediction Score")
    print(p1)
  } else{
    p2 <- ggplot(cell_meta, aes(pred_score, fill = annotation_correct, colour = annotation_correct)) +
      geom_density(position = "stack",aes(y = after_stat(count),alpha = 0.5)) +
      scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      xlab("Prediction Score")
    print(p2)
  }
}










plot_central_cells <- function(Seurat.object, cor_mtx, label = T, pt.size = 5,reduction = DefaultDimReduc(Seurat.object)){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()

  central_cell_idx <- apply(cor_mtx,2,which.max)
  Seurat.object$test <- ' '
  Seurat.object$test[central_cell_idx] <- names(central_cell_idx)
  p <- plot_highlight_cells(Seurat.object,setdiff(unique(Seurat.object$test),' '),ident = 'test',
                            label = label,pt.size = pt.size,reduction = reduction) + NoLegend()
  return(p)
}





