################### ResultsPlotFunctions ##################
# plot_confusion_matrix
# plot_highlight_cells
# plot_cell_distribution
# plot_pred_scores
# plot_seed_cells

#' Plot highlight cells
#'
#' @param Seurat.object A Seurat object
#' @param celltype Which cell type(s) to highlight
#' @param ident The identity class to plot
#' @param label Whether to display the cell type label on the plot, default is TRUE
#' @param pt.size Point size, default is 1
#' @param reduction Which dimensionality reduction to use
#' @param colors Colors of different cell types
#'
#' @return Return a ggplot object
#' @export
#'
plot_highlight_cells <- function(Seurat.object, celltype, ident, label = T, pt.size = 1,
                                 reduction = DefaultDimReduc(Seurat.object),colors = NULL){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  Idents(Seurat.object) <- ident
  cells.highlight <- list()
  for (i in celltype){
    cells.highlight[[`i`]] <- WhichCells(Seurat.object, idents = i)
  }
  if(is.null(colors)){
    colors <- c("red", "orange",'yellow','green','skyblue','pink2','purple','darkgreen','blue')
  }
  p <- DimPlot(Seurat.object, label = label, group.by = ident, cells.highlight = cells.highlight, reduction = reduction,
               cols.highlight = colors,
               sizes.highlight = pt.size,repel = T) + ggtitle(ident)
  return(p)
}




#' Plot global markers heatmap
#'
#' @param ref_mtx The reference gene expression matrix
#' @param ref_labels The cell type labels of the reference cells
#' @param global_markers `global_markers` generate from the function `get_global_markers_()`
#' @param neighbor_celltypes `neighbor_celltypes` generate from the function `get_neighbor_celltypes()`
#' @param celltypes_to_plot Cell type to display, default is all cell types
#' @param top_marker_genes_num Number of top marker genes to display for each cell type, default is 20
#' @param sample_cells_num Number of the sampling cells to display for each cell type, default is 20
#'
#' @export
#'
plot_ref_global_markers_heatmap <- function(ref_mtx,ref_labels,
                                            global_markers,neighbor_celltypes = NULL,celltypes_to_plot = NULL,
                                            top_marker_genes_num = 20,sample_cells_num = 20){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  if(is.null(celltypes_to_plot)){
    cts <- names(neighbor_celltypes)
  } else{
    cts <- celltypes_to_plot
  }
  ordered_genes <- sapply(global_markers[cts],function(x){x$global_markers[1:min(length(x$global_markers),top_marker_genes_num)]}) %>% 
    unlist() %>% as.vector()
  all_ct_idx <- c()
  for(ct in cts){
    # ct <- 'Plasma'
    ct_idx <- which(ref_labels == ct)
    sample_num <- min(length(ct_idx),sample_cells_num)
    ct_sample_idx <- sample(ct_idx,sample_num)
    all_ct_idx <- c(all_ct_idx,ct_sample_idx)
  }
  ref_mtx_norm <- Seurat::NormalizeData(ref_mtx,verbose = F)
  heat_data <- ref_mtx_norm[ordered_genes,all_ct_idx] %>% as.matrix()
  Heatmap(t(scale(t(scale(heat_data)))), cluster_rows = F, cluster_columns = F, 
          name = 'Scaled\nexpression',use_raster = F,
          row_title = "Global marker genes",column_title = "Cells",
          show_column_names = F, show_row_names = F, 
          bottom_annotation = HeatmapAnnotation(celltype = factor(ref_labels[all_ct_idx],levels = cts))) %>% 
    suppressMessages %>% print
}





#' Plot neighbor markers heatmap
#'
#' @param ref_mtx The reference gene expression matrix
#' @param ref_labels The cell type labels of the reference cells
#' @param neighbor_markers `neighbor_markers` generate from the function `get_neighbor_markers_()`
#' @param neighbor_celltypes `neighbor_celltypes` generate from the function `get_neighbor_celltypes()`
#' @param celltypes_to_plot Which cell type and its neighbor cell types to display
#' @param top_marker_genes_num Number of top marker genes to display for each cell type, default is 20
#' @param sample_cells_num Number of the sampling cells to display for each cell type, default is 20
#'
#' @export
#'
plot_ref_neighbor_markers_heatmap <- function(ref_mtx,ref_labels,
                                              neighbor_markers,neighbor_celltypes = NULL,celltype_to_plot,
                                              top_marker_genes_num = 20,sample_cells_num = 20){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  cts <- neighbor_celltypes[[celltype_to_plot]]
  neighbor_marker_list <- lapply(neighbor_markers[cts],function(x){x$neighbor_markers[1:min(length(x$neighbor_markers),top_marker_genes_num)]}) 
  ordered_genes <- neighbor_marker_list %>% unlist() %>% as.vector()
  all_ct_idx <- c()
  for(ct in cts){
    # ct <- 'Plasma'
    # print(ct)
    ct_idx <- which(ref_labels == ct)
    sample_num <- min(length(ct_idx),sample_cells_num)
    ct_sample_idx <- sample(ct_idx,sample_num)
    all_ct_idx <- c(all_ct_idx,ct_sample_idx)
  }
  
  ref_mtx_norm <- Seurat::NormalizeData(ref_mtx,verbose = F)
  heat_data <- ref_mtx_norm[ordered_genes,all_ct_idx] %>% as.matrix()
  data <- t(scale(t(heat_data)))
  Heatmap(data, cluster_rows = F, cluster_columns = F, 
          show_column_names = F, show_row_names = F, 
          row_title = "Neighbor marker genes",use_raster = F,
          name = 'Scaled\nexpression',
          row_split = rep(names(neighbor_marker_list),sapply(neighbor_marker_list,length)), 
          column_split = ref_labels[all_ct_idx],
          bottom_annotation = HeatmapAnnotation(celltype = factor(ref_labels[all_ct_idx],levels = cts)))  %>% 
    suppressMessages %>% print
}






#' Plot global markers heatmap of seed cells gene activity
#'
#' @param query_mtx The query gene activity matrix
#' @param cell_meta The cell metadata
#' @param global_markers `global_markers` generate from the function `get_global_markers_()`
#' @param neighbor_celltypes `neighbor_celltypes` generate from the function `get_neighbor_celltypes()`
#' @param celltypes_to_plot Cell type to display, default is all cell types
#' @param top_marker_genes_num Number of top marker genes to display for each cell type, default is 20
#' @param sample_cells_num Number of the sampling cells to display for each cell type, default is 20
#'
#' @return Return a Heatmap-class object
#' @export
#'
plot_seed_global_markers_heatmap <- function(query_mtx,cell_meta,
                                             global_markers,neighbor_celltypes,celltypes_to_plot = NULL,
                                             top_marker_genes_num = 20,sample_cells_num = 20){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  if(is.null(celltypes_to_plot)){
    cts <- names(neighbor_celltypes)
  } else{
    cts <- celltypes_to_plot
  }
  ordered_genes <- sapply(global_markers[cts],function(x){x$global_markers[1:min(length(x$global_markers),top_marker_genes_num)]}) %>% 
    unlist() %>% as.vector()
  all_ct_idx <- c()
  for(ct in cts){
    # ct <- 'Plasma'
    ct_idx <- which(seed_meta$kendall_pred == ct)
    sample_num <- min(length(ct_idx),20)
    ct_sample_idx <- sample(ct_idx,sample_num)
    all_ct_idx <- c(all_ct_idx,ct_sample_idx)
  }
  all_ct_barcodes <- rownames(seed_meta)[all_ct_idx]
  query_mtx_norm <- Seurat::NormalizeData(query_mtx,verbose = F)
  heat_data <- query_mtx_norm[ordered_genes,all_ct_barcodes] %>% as.matrix()
  data <- t(scale(t(heat_data)))
  Heatmap(data, cluster_rows = F, cluster_columns = F, use_raster = F,
          row_title = "Global marker genes",column_title = "Cells",
          name = 'Scaled\ngene activity',
          show_column_names = F, show_row_names = F, 
          bottom_annotation = HeatmapAnnotation(celltype = factor(seed_meta$kendall_pred[all_ct_idx],levels = cts))) %>% 
    suppressMessages %>% print()
}





#' Plot neighbor markers heatmap of seed cells gene activity
#'
#' @param query_mtx The query gene activity matrix
#' @param cell_meta The cell metadata
#' @param neighbor_markers `neighbor_markers` generate from the function `get_neighbor_markers_()`
#' @param neighbor_celltypes `neighbor_celltypes` generate from the function `get_neighbor_celltypes()`
#' @param celltypes_to_plot Cell type to display, default is all cell types
#' @param top_marker_genes_num Number of top marker genes to display for each cell type, default is 20
#' @param sample_cells_num Number of the sampling cells to display for each cell type, default is 20
#'
#' @return Return a Heatmap-class object
#' @export
#'
plot_seed_neighbor_markers_heatmap <- function(query_mtx,cell_meta,
                                               neighbor_markers = NULL,neighbor_celltypes = NULL,celltype_to_plot,
                                               top_marker_genes_num = 20,sample_cells_num = 20){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  cts <- neighbor_celltypes[[celltype_to_plot]]
  neighbor_marker_list <- lapply(neighbor_markers[cts],function(x){x$neighbor_markers[1:min(length(x$neighbor_markers),top_marker_genes_num)]}) 
  ordered_genes <- neighbor_marker_list %>% unlist() %>% as.vector()
  all_ct_idx <- c()
  for(ct in cts){
    # ct <- 'Plasma'
    ct_idx <- which(seed_meta$kendall_pred == ct)
    sample_num <- min(length(ct_idx),20)
    ct_sample_idx <- sample(ct_idx,sample_num)
    all_ct_idx <- c(all_ct_idx,ct_sample_idx)
  }
  all_ct_barcodes <- rownames(seed_meta)[all_ct_idx]
  
  heat_data <- Seurat::NormalizeData(query_mtx[ordered_genes,all_ct_barcodes],verbose = F) %>% as.matrix()
  data <- t(scale(t(heat_data)))
  Heatmap(data, cluster_rows = F, cluster_columns = F, use_raster = F,
          show_column_names = F, show_row_names = F, 
          name = 'Scaled\ngene activity',row_title = "Neighbor marker genes",
          row_split = rep(names(neighbor_marker_list),sapply(neighbor_marker_list,length)), 
          column_split = seed_meta$kendall_pred[all_ct_idx],
          bottom_annotation = HeatmapAnnotation(celltype = factor(seed_meta$kendall_pred[all_ct_idx],levels = cts))) %>% 
    suppressMessages %>% print()
}








#' Plot the UMAP of the predictions
#'
#' @param Seurat.object A Seurat Object
#' @param cell_meta The cell metadata
#' @param celltype which cell types to show
#' @param category which predictions to show, either `seed_candidate`, `seed`, or `final`, default is `final`
#' @param label Whether to display the cell type label on the plot, default is TRUE
#' @param pt.size Point size, default is 1
#' @param reduction Which dimensionality reduction to use
#' @param colors Colors of different cell types
#'
#' @return Return a ggplot object
#' @export
#'
plot_pred_umap <- function(Seurat.object, cell_meta, celltype = NULL,category = 'final',
                           label = T,pt.size = 0.3,reduction = DefaultDimReduc(Seurat.object),colors = NULL){
  library(ggplot2) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  library(Seurat) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  if(category == 'seed'){
    final_seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
    final_seed_barcodes <- rownames(final_seed_meta)
    
    Seurat.object$seed <- ' '
    Seurat.object@meta.data[final_seed_barcodes,]$seed <- final_seed_meta$kendall_pred
    if(is.null(celltype)){
      p <- plot_highlight_cells(Seurat.object,celltype = setdiff(unique(Seurat.object$seed),' '),'seed',
                                pt.size = pt.size,reduction = reduction,colors = colors,label = label) + 
        ggtitle('Seed cells')
    } else{
      p <- plot_highlight_cells(Seurat.object,celltype = celltype,'seed',
                                pt.size = pt.size,reduction = reduction,colors = colors,label = label) + 
        ggtitle('Seed cells')
    }
  } else if(category == 'seed_candidate'){
    final_seed_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
    final_seed_barcodes <- rownames(final_seed_meta)
    
    Seurat.object$seed_candidate <- ' '
    Seurat.object@meta.data[final_seed_barcodes,]$seed_candidate <- final_seed_meta$kendall_pred
    if(is.null(celltype)){
      p <- plot_highlight_cells(Seurat.object,celltype = setdiff(unique(Seurat.object$seed_candidate),' '),
                                'seed_candidate',pt.size = pt.size,colors = colors,label = label) + 
        ggtitle('Seed cell candidates')
    } else{
      p <- plot_highlight_cells(Seurat.object,celltype = celltype,'seed_candidate',
                                label = label,pt.size = pt.size,colors = colors) + 
        ggtitle('Seed cell candidates')
    } 
  } else{
    Seurat.object$pred <- cell_meta$final_pred
    if(is.null(celltype)){
      p <- DimPlot(Seurat.object,group.by = 'pred',label = label,pt.size = pt.size,cols = colors) + 
        ggtitle('Final predictions')
    } else{
      p <- plot_highlight_cells(Seurat.object,celltype = celltype,ident = 'pred',label = label,pt.size = pt.size,colors = colors) + 
        ggtitle('Final predictions')
    }
  }
  return(p)
}





good_heatmap <- function(data,label){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  # label <- seed_meta$kendall_pred[all_cell_idx]
  cts <- unique(label)
  nmf_pb <- get_pseudo_bulk_mtx(data,labels = label)[,cts]
  ideal_markers <- sapply(1:ncol(nmf_pb),function(row_idx){
    row <- nmf_pb[row_idx,]
    ideal_marker <- rep(0,ncol(nmf_pb))
    ideal_marker[row_idx] <- 1
    ideal_marker
  })
  cos_mtx <- RcppML::cosine(t(nmf_pb),ideal_markers)
  dimnames(cos_mtx) <- list(rownames(data),colnames(nmf_pb))
  factor_1to50_cts <- colnames(nmf_pb)[apply(cos_mtx,1,which.max)]
  names(factor_1to50_cts) <- rownames(cos_mtx)
  ordered_factors <- sort(factor(factor_1to50_cts,levels = cts)) %>% names()
  
  Heatmap(data[ordered_factors,], cluster_rows = F, cluster_columns = F, name = 'Cell\nembedding',
          row_title = "Meta-programs",column_title = "Cells",
          show_column_names = F, show_row_names = T, use_raster = F,row_names_gp = gpar(fontsize = 8),
          bottom_annotation = HeatmapAnnotation(celltype = factor(label,levels = cts))) %>% 
    suppressMessages %>% return
}







#' Plot the NMF meta-program heatmap of the predictions
#'
#' @param query_nmf_embedding The query nmf embedding
#' @param cell_meta The cell metadata
#' @param neighbor_celltypes `neighbor_celltypes` generate from the function `get_neighbor_celltypes()`
#' @param category which predictions to show, either `'seed_candidate'`, `'seed'`, or `final`, default is `'final'`
#' @param sample_cells Whether to sample cells, default is FALSE, if is true, randomly sample 30 cells for each cell type
#' @param celltypes_to_plot Cell type to display, default is all cell types
#'
#' @return Return a Heatmap-class object
#' @export
#'
plot_pred_nmf <- function(query_nmf_embedding,cell_meta,neighbor_celltypes,category = 'final',
                          sample_cells = F,celltypes_to_plot = NULL){
  library(ComplexHeatmap) %>% suppressPackageStartupMessages() %>% suppressMessages() %>% suppressWarnings()
  
  if(is.null(celltypes_to_plot)){
    cts <- names(neighbor_celltypes)
  } else{
    cts <- celltypes_to_plot
  }
  
  if(category == 'seed'){
    seed_idx <- which(cell_meta$is_seed == T)
  } else if(category == 'seed_candidate'){
    seed_idx <- which(cell_meta$is_seed_candidate == T)
  } else{
    seed_idx <- 1:nrow(cell_meta)
  }
  
  seed_nmf <- query_nmf_embedding[seed_idx,]
  seed_meta <- cell_meta[seed_idx,]
  
  all_cell_idx <- c()
  for (ct in cts){
    ct_cell_idx <- which(seed_meta$kendall_pred == ct)
    if(sample_cells){
      sample_ct_cell_idx <- sample(ct_cell_idx,min(30,length(ct_cell_idx)))
      all_cell_idx <- c(all_cell_idx,sample_ct_cell_idx)
    } else{
      all_cell_idx <- c(all_cell_idx,ct_cell_idx)
    }
  }
  ordered_seed_nmf <- seed_nmf[all_cell_idx,]
  data <- t(ordered_seed_nmf)
  good_heatmap(data,label = seed_meta$kendall_pred[all_cell_idx]) %>% print()
}






#' Plot cell type proportions for the reference, candidate seed cells, seed cells and the final predictions
#'
#' @param cell_meta The cell metadata
#' @param ref_labels The cell type labels of the reference cells
#'
#' @return Return a ggplot object
#' @export
#'
plot_celltype_proportions <- function(cell_meta,ref_labels){
  
  reference_df <- as.data.frame(table(unname(ref_labels)))
  reference_df$class <- 'reference'
  
  seed_candidate_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
  seed_candidate_meta_df <- as.data.frame(table(seed_candidate_meta$kendall_pred))
  seed_candidate_meta_df$class <- 'seed_cell_candidates'
  
  seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  seed_meta_df <- as.data.frame(table(seed_meta$kendall_pred))
  seed_meta_df$class <- 'seed_cells'
  
  cell_meta_df <- as.data.frame(table(cell_meta$final_pred))
  cell_meta_df$class <- 'final_pred'
  
  plot_df <- rbind(reference_df,seed_candidate_meta_df,seed_meta_df,cell_meta_df)
  plot_df$class <- factor(plot_df$class,
                          levels = c('reference','seed_cell_candidates','seed_cells','final_pred') %>% rev(),
                          labels = c('Reference','Seed\ncell\ncandidates','Seed\ncells','Final\npredictions') %>% rev())
  ggplot(plot_df,mapping = aes(Freq,class,fill=Var1))+
    geom_bar(stat='identity',position='fill',color = 'black') +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14, color = 'black'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(color = 'black',size = 10)) + 
    theme(axis.text.y = element_text(color = 'black',size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    # theme(panel.border = element_blank()) +
    labs(x = 'Proportion',y = '',fill = 'Cell types') +
    scale_y_discrete(position = "right") 
  
}







