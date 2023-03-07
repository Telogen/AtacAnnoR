################### PlotFunctions ##################
# plot_confusion_matrix
# plot_highlight_cells
# plot_cell_distribution
# plot_pred_scores
# plot_seed_cells

plot_confusion_matrix <- function(cell_meta, which.group, mode = 1,title = 3){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(cowplot) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
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
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
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
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
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
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # +
    #   geom_text(aes(label = Freq),size = 3)
    
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
                                   x = 'NMSS',y = 'GMSS',print = T){
  library(aplot) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  if(is.null(cell_meta[['is_seed']])){
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
  P_list <- list()
  for(i in cell_meta_li){
    DATA <- i
    
    p1 <- ggplot(data = DATA) +
      geom_point(aes_string(x = x,y = y, color = name_booltrue,shape = 'is_seed')) +
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







plot_seed_cells <- function(Seurat.object, cell_meta, celltype = NULL,reduction = DefaultDimReduc(Seurat.object),pt.size = 0.3){
  final_seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  final_seed_barcodes <- rownames(final_seed_meta)

  Seurat.object$seed <- ' '
  Seurat.object@meta.data[final_seed_barcodes,]$seed <- final_seed_meta$kendall_pred
  if(is.null(celltype)){
    p <- plot_highlight_cells(Seurat.object,celltype = setdiff(unique(Seurat.object$seed),' '),'seed',pt.size = pt.size) + NoLegend()
  } else{
    p <- plot_highlight_cells(Seurat.object,celltype = celltype,'seed',pt.size = pt.size) + NoLegend()
  }
  return(p)
}


