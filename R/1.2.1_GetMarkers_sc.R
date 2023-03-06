########### GetMarkers_sc ###########
# get_global_markers_sc
# get_neighbor_markers_sc




#' Get global markers for scRNA-seq gene counts matrix
#' 
#' Get global markers for scRNA-seq gene counts matrix
#'
#' @param sc_counts_mtx scRNA-seq counts matrix
#' @param labels scRNA-seq cell labels
#' @param max_marker maximum number of markers for each cell type
#' @param threads the number of threads, default is 10
#' @param return_raw whether return raw data, default is FALSE
#'
#' @return Returns a list of global markers for each cell type
#' @export
get_global_markers_sc <- function(sc_counts_mtx, labels,max_marker = 200,threads = 10,return_raw = F){
  # sc_counts_mtx <- ref_mtx
  # labels <- SeuratObj_RNA$true
  # return.raw = F
  
  # if one ct more than 2000 cells, sample 2000 cells
  if(max(table(labels)) > 2000){
    label_list <- split(1:ncol(sc_counts_mtx),labels)
    sampled_idx <- lapply(label_list,function(label_i){
      if(length(label_i) > 2000){
        return(sample(label_i,2000))
      } else{
        return(label_i)
      }
    }) %>% 
      unlist(use.names = F) %>%
      sort()
    sc_counts_mtx <- sc_counts_mtx[,sampled_idx]
    labels <- labels[sampled_idx]
  }
  
  # creat Seurat object
  sc_counts_mtx_Seurat <- CreateSeuratObject(sc_counts_mtx) %>% NormalizeData(verbose = F)
  sc_counts_mtx_Seurat$true <- labels
  all_cts <- names(table(sc_counts_mtx_Seurat$true))
  Idents(sc_counts_mtx_Seurat) <- 'true'
  
  # get each ct's highly expressed genes
  pb_ref <- get_pseudo_bulk_mtx(sc_counts_mtx, labels)
  HEG_num <- round(nrow(pb_ref) * 0.2)
  HEG_df <- apply(pb_ref,2,function(col){rownames(pb_ref)[order(col,decreasing = T)][1:HEG_num]}) %>% as.data.frame()
  
  # get each ct's markers
  Seurat_marker_list <- pbmcapply::pbmclapply(all_cts,function(ct){
    HEG <- HEG_df[[ct]]
    mak_df <- Seurat::FindMarkers(sc_counts_mtx_Seurat,features = HEG,ident.1 = ct,
                                  only.pos = T,logfc.threshold = 0.5,min.pct = 0,min.cells.group = 0,verbose = F)
    mak_df$p_val_adj <- stats::p.adjust(p = mak_df$p_val,method = 'fdr')
    mak_df$pct_fc <- mak_df$pct.1/mak_df$pct.2
    mak_df$pct_diff <- mak_df$pct.1 - mak_df$pct.2
    mak_df <- dplyr::filter(mak_df,p_val_adj < 0.01, pct_fc > 1.5 | pct_diff > 0.1) %>% 
      dplyr::arrange(desc(avg_log2FC))
    if(nrow(mak_df) > max_marker){
      markers <- rownames(mak_df)[1:max_marker]
    } else{
      markers <- rownames(mak_df)
    }
    if(return_raw){
      return(mak_df)
    } else{
      return(markers)
    }
  },mc.cores = threads)
  names(Seurat_marker_list) <- all_cts
  
  
  if(return_raw){
    return(Seurat_marker_list)
  } else{
    # get global_bg_genes for each ct
    global_markers <- list()
    for(i in 1:length(Seurat_marker_list)){
      i_markers <- Seurat_marker_list[[i]]
      other_cts_markers <- Seurat_marker_list[-i] %>% unlist(use.names = F) %>% unique()
      global_markers[[i]] <- list(global_markers = i_markers,
                                  global_bg_genes = setdiff(other_cts_markers,i_markers))
    }
    names(global_markers) <- all_cts
    return(global_markers)
  }
}






#' Get neighbor markers for scRNA-seq gene counts matrix
#' 
#' Get neighbor markers for scRNA-seq gene counts matrix
#'
#' @param sc_counts_mtx scRNA-seq counts matrix
#' @param labels scRNA-seq cell labels
#' @param neighbor_celltypes neighbor celltypes got from `get_neighbor_celltypes()`
#' @param global_markers global markers got from `get_global_markers_sc()`
#' @param max_marker maximum number of markers for each cell type
#' @param threads the number of threads, default is 10
#'
#' @return Returns a list of neighbor markers for each cell type
#' @export
get_neighbor_markers_sc <- function(sc_counts_mtx, labels, neighbor_celltypes, global_markers,max_marker = 200,threads = 10){
  # sc_counts_mtx <- ref_mtx
  # labels <- SeuratObj_RNA$true
  
  # if one ct more than 2000 cells, sample 2000 cells
  if(max(table(labels)) > 2000){
    label_list <- split(1:ncol(sc_counts_mtx),labels)
    sampled_idx <- lapply(label_list,function(label_i){
      if(length(label_i) > 2000){
        return(sample(label_i,2000))
      } else{
        return(label_i)
      }
    }) %>% 
      unlist(use.names = F) %>%
      sort()
    sc_counts_mtx <- sc_counts_mtx[,sampled_idx]
    labels <- labels[sampled_idx]
  }
  
  # creat Seurat object
  sc_counts_mtx_Seurat <- CreateSeuratObject(sc_counts_mtx) %>% NormalizeData(verbose = F)
  sc_counts_mtx_Seurat$true <- labels
  all_cts <- names(table(sc_counts_mtx_Seurat$true))
  Idents(sc_counts_mtx_Seurat) <- 'true'
  
  # get each ct's highly expressed genes
  pb_ref <- get_pseudo_bulk_mtx(sc_counts_mtx, labels)
  HEG_num <- round(nrow(pb_ref) * 0.2)
  HEG_df <- apply(pb_ref,2,function(col){rownames(pb_ref)[order(col,decreasing = T)][1:HEG_num]}) %>% as.data.frame()
  
  # get each ct's neighbor markers
  Seurat_marker_list <- pbmcapply::pbmclapply(all_cts,function(ct){
    neighbor_cts <- neighbor_celltypes[[ct]][-1]
    
    if(length(neighbor_cts) == 0){
      markers <- NA
    } else{
      HEG <- HEG_df[[ct]]
      mak_df <- Seurat::FindMarkers(sc_counts_mtx_Seurat,features = HEG,ident.1 = ct,ident.2 = neighbor_cts,
                                    only.pos = T,logfc.threshold = 0.5,min.pct = 0,min.cells.group = 0,verbose = F)
      mak_df$p_val_adj <- stats::p.adjust(p = mak_df$p_val,method = 'fdr')
      mak_df$pct_fc <- mak_df$pct.1/mak_df$pct.2
      mak_df$pct_diff <- mak_df$pct.1 - mak_df$pct.2
      mak_df <- dplyr::filter(mak_df,p_val_adj < 0.01, pct_fc > 1.5 | pct_diff > 0.1) %>% 
        dplyr::arrange(desc(avg_log2FC))
      
      if(nrow(mak_df) < 20){
        mak_df <- Seurat::FindMarkers(sc_counts_mtx_Seurat,features = HEG,ident.1 = ct,ident.2 = neighbor_cts,
                                      only.pos = T,logfc.threshold = 0.1,min.pct = 0,min.cells.group = 0,verbose = F)
        mak_df$p_val_adj <- stats::p.adjust(p = mak_df$p_val,method = 'fdr')
        mak_df <- dplyr::filter(mak_df,p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC))
      }
      
      if(nrow(mak_df) > max_marker){
        markers <- rownames(mak_df)[1:max_marker]
      } else{
        markers <- rownames(mak_df)
      }
    }
    return(markers)
  },mc.cores = threads)
  names(Seurat_marker_list) <- all_cts
  
  # get neighbor_bg_genes for each ct
  neighbor_markers <- list()
  for(ct in all_cts){
    if(length(neighbor_celltypes[[ct]]) == 1){
      neighbor_markers[[ct]]$neighbor_markers <- global_markers[[ct]]$global_markers
      neighbor_markers[[ct]]$neighbor_bg_genes <- global_markers[[ct]]$global_bg_genes
    } else{
      neighbor_cts <- neighbor_celltypes[[ct]][-1]
      neighbor_bg_genes_tmp <- unlist(Seurat_marker_list[neighbor_cts],use.names = F) %>% unique()
      neighbor_markers[[ct]]$neighbor_markers <- Seurat_marker_list[[ct]]
      neighbor_markers[[ct]]$neighbor_bg_genes <- setdiff(neighbor_bg_genes_tmp,Seurat_marker_list[[ct]])
    }
  }
  
  return(neighbor_markers)
}







