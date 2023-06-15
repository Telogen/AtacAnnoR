################### BenchmarkPlotFunctions ##################
# plot_confusion_matrix
# plot_highlight_cells
# plot_cell_distribution
# plot_pred_scores
# plot_seed_cells


plot_confusion_matrix <- function(cell_meta, which.group, mode = 1,title = 3){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()

  bench <- get_benchmark(cell_meta$true,cell_meta[,which(colnames(cell_meta) == which.group)])
  bench <- round(bench*100,2)
  # if(is.null(title)){
  #   title <- mode
  # }
  
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
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,color = 'black')) +
      theme(axis.text.y = element_text(color = 'black')) + 
      theme(title = element_text(hjust = 0.5,color  = 'black'))
    # +
    #   geom_text(aes(label = freq_morethan5),size = 3)
    
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
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,color = 'black')) +
      theme(axis.text.y = element_text(color = 'black')) + 
      theme(title = element_text(hjust = 0.5,color  = 'black'))
    # +
    #   geom_text(aes(label = freq_morethan5),size = 3)
    
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
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,color = 'black')) +
      theme(axis.text.y = element_text(color = 'black')) + 
      theme(title = element_text(hjust = 0.5,color  = 'black'))
    # +
    #   geom_text(aes(label = Freq),size = 3)
    
  }
  
  return(g)
}



plot_cell_distribution <- function(cell_meta,celltypes.to.plot = NULL,group.by = 'kendall_pred',
                                   x = 'NMSS',y = 'GMSS',pt.size = 2,print = T){
  library(aplot) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  if(is.null(cell_meta[['is_seed']])){
    cell_meta$is_seed <- FALSE
  }
  name_booltrue <- 'Correctness'
  cell_meta[[name_booltrue]] <- factor((cell_meta[,group.by] == cell_meta[,'true']))
  
  cell_meta_li <- split(cell_meta,cell_meta[,group.by])
  
  if (is.null(celltypes.to.plot)){
    cell_meta_li <- cell_meta_li
  } else{
    cell_meta_li <- cell_meta_li[celltypes.to.plot]
  }
  P_list <- list()
  for(i in cell_meta_li){
    DATA <- i
    
    p1 <- ggplot(data = DATA) +
      geom_point(aes_string(x = x,y = y, color = name_booltrue,shape = 'is_seed'),size = pt.size) +
      scale_x_continuous(expand=c(0,0)) + 
      scale_y_continuous(expand=c(0,0))
    
    
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
    
    if(!is.null(cell_meta$GMSS_cutoff)){
      p1 <- p1 + geom_hline(yintercept = DATA$GMSS_cutoff[1], linetype = 2)
    }
    if(!is.null(cell_meta$NMSS_cutoff)){
      seed_candidate_DATA <- DATA[which(DATA$is_seed_candidate),]
      seed_DATA <- DATA[which(DATA$is_seed),]
      p1 <- p1 + 
        geom_vline(xintercept = DATA$NMSS_cutoff[1], linetype = 2)
      # + 
      #   annotate("text",x=Inf, y = Inf, vjust=2, hjust=1,size = 5,
      #            label = paste0(' Kept cell number: ',
      #                           nrow(seed_candidate_DATA),
      #                           ' (',
      #                           round(nrow(seed_candidate_DATA)/nrow(DATA)*100,2),
      #                           '%) \n Seed candidate ccuracy: ',
      #                           round(get_benchmark(seed_candidate_DATA$true,seed_candidate_DATA$kendall_pred)[1],3),
      #                           ' \n Final seed number: ',
      #                           nrow(seed_DATA),
      #                           ' (',
      #                           round(nrow(seed_DATA)/nrow(DATA)*100,2),
      #                           '%) \n Final seed accuracy: ',
      #                           round(get_benchmark(seed_DATA$true,seed_DATA$kendall_pred)[1],3)
      #                           )
      #            )
    }
    
    P <- p1 %>% insert_top(p2,height = 0.2) %>% insert_right(p3,0.2)
    
    if (print){
      print(P)
    } else{
      ct <- DATA$kendall_pred[1]
      P_list[[ct]] <- P
    }
  }
  
  if(!print){
    return(P_list)
  }
}






plot_cell_correctness <- function(Seurat.object,cell_meta,which_cells,
                                  pt.size = 0.5,label = F,reduction = DefaultDimReduc(Seurat.object)){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  cell_meta$kendall_pred_booltrue <- cell_meta$kendall_pred == cell_meta$true
  if(which_cells == 'seed'){
    seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
    seed_barcodes <- rownames(seed_meta)
    Seurat.object$cell_correctness <- ' '
    Seurat.object@meta.data[seed_barcodes,]$cell_correctness <- seed_meta$kendall_pred_booltrue
    cells.highlight <- list('TRUE' = which(Seurat.object$cell_correctness == T),
                            'FALSE' = which(Seurat.object$cell_correctness == F))
    
    p <- DimPlot(Seurat.object, label = label, group.by = 'cell_correctness', 
                 cells.highlight = cells.highlight, reduction = reduction,
                 cols.highlight = c('red','blue'),
                 sizes.highlight = pt.size,repel = T) + 
      ggtitle(paste0('Seed cells accuracy: ',
                     round(sum(seed_meta$kendall_pred_booltrue)/nrow(seed_meta)*100,2),'%'))  
  } else{
    seed_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
    seed_barcodes <- rownames(seed_meta)
    Seurat.object$cell_correctness <- ' '
    Seurat.object@meta.data[seed_barcodes,]$cell_correctness <- seed_meta$kendall_pred_booltrue
    cells.highlight <- list('TRUE' = which(Seurat.object$cell_correctness == T),
                            'FALSE' = which(Seurat.object$cell_correctness == F))
    
    p <- DimPlot(Seurat.object, label = label, group.by = 'cell_correctness', 
                 cells.highlight = cells.highlight, reduction = reduction,
                 cols.highlight = c('red','blue'),
                 sizes.highlight = pt.size,repel = T) + 
      ggtitle(paste0('Seed cell candidates accuracy: ',
                     round(sum(seed_meta$kendall_pred_booltrue)/nrow(seed_meta)*100,2),'%'))  
  }
  return(p)
}





