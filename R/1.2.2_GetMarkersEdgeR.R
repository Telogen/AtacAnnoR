########### GetMarkersEdgeR ###########
# get_global_markers_edgeR
# get_neighbor_markers_edgeR





get_global_markers_edgeR <- function(sc_counts_mtx,labels,max_marker = 200,threads = 10,return_raw = F){
  
  pb_ref <- get_pseudo_bulk_mtx(sc_counts_mtx, labels)
  HEG_num <- round(nrow(pb_ref) * 0.5)
  HEG_df <- apply(pb_ref,2,function(col){rownames(pb_ref)[order(col,decreasing = T)][1:HEG_num]}) %>% as.data.frame()
  
  all_cts <- names(table(labels))
  
  edgeR_marker_list <- pbmcapply::pbmclapply(all_cts,function(ct){
    HEG <- HEG_df[[ct]]
    group <- factor(as.numeric(labels == ct))
    design <- model.matrix(~ 0 + group)
    dge <- edgeR::DGEList(sc_counts_mtx, group=group)
    keep.exprs <- edgeR::filterByExpr(dge) #自动筛选过滤低表达基因
    dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
    dge <- edgeR::calcNormFactors(dge, method = 'TMM') #归一化因子用于 normalizes the library sizes
    dge <- edgeR::estimateDisp(dge, design, robust=T) 
    fit <- edgeR::glmFit(dge, design, robust=T)
    lt <- edgeR::glmLRT(fit, contrast=c(-1,1))  #比对:顺序实验/对照，已设对照为1
    tempDEG <- edgeR::topTags(lt, n = Inf,sort.by = 'logFC') 
    tempDEG <- as.data.frame(tempDEG)
    DEG_edgeR <- na.omit(tempDEG)
    DEG_edgeR$gene <- rownames(DEG_edgeR)
    
    select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.5,FDR < 0.05,gene %in% HEG)
    nrow(select_DEG_edgeR)
    if(nrow(select_DEG_edgeR) < 50){
      select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.2,PValue < 0.01)
    }
    
    if(nrow(select_DEG_edgeR) > 200){
      select_DEG_edgeR <- select_DEG_edgeR[1:200,]
    }
    
    markers_edgeR <- rownames(select_DEG_edgeR)
    
    if(return_raw){
      colnames(select_DEG_edgeR)[1] <- 'avg_log2FC'
      return(select_DEG_edgeR)
    } else{
      return(markers_edgeR)
    }

  },mc.cores = threads)
  names(edgeR_marker_list) <- all_cts
  # sapply(edgeR_marker_list,length)
  
  if(return_raw){
    
    return(edgeR_marker_list)
  } else{
    # get global_bg_genes for each ct
    global_markers <- list()
    for(i in 1:length(edgeR_marker_list)){
      i_markers <- edgeR_marker_list[[i]]
      other_cts_markers <- edgeR_marker_list[-i] %>% unlist(use.names = F) %>% unique()
      global_markers[[i]] <- list(global_markers = i_markers,
                                  global_bg_genes = setdiff(other_cts_markers,i_markers))
    }
    names(global_markers) <- all_cts
    return(global_markers)
  }
  
}






get_neighbor_markers_edgeR <- function(sc_counts_mtx,labels,neighbor_celltypes, global_markers,max_marker = 200,threads = 10){
  
  pb_ref <- get_pseudo_bulk_mtx(sc_counts_mtx, labels)
  HEG_num <- round(nrow(pb_ref) * 0.5)
  HEG_df <- apply(pb_ref,2,function(col){rownames(pb_ref)[order(col,decreasing = T)][1:HEG_num]}) %>% as.data.frame()
  
  all_cts <- names(table(labels))
  
  edgeR_marker_list <- pbmcapply::pbmclapply(all_cts,function(ct){
    neighbor_cts <- neighbor_celltypes[[ct]]
    
    if(length(neighbor_cts) == 1){
      markers_edgeR <- NA
    } else{
      HEG <- HEG_df[[ct]]
      ct_neighbor_idx <- which(labels %in% neighbor_cts)
      labels2 <- labels[ct_neighbor_idx]
      sc_counts_mtx2 <- sc_counts_mtx[,ct_neighbor_idx]
      
      group <- factor(as.numeric(labels2 == ct))
      design <- model.matrix(~ 0 + group)
      dge <- edgeR::DGEList(sc_counts_mtx2, group=group)
      keep.exprs <- edgeR::filterByExpr(dge) #自动筛选过滤低表达基因
      dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
      dge <- edgeR::calcNormFactors(dge, method = 'TMM') #归一化因子用于 normalizes the library sizes
      dge <- edgeR::estimateDisp(dge, design, robust=T) 
      fit <- edgeR::glmFit(dge, design, robust=T)
      lt <- edgeR::glmLRT(fit, contrast=c(-1,1))  #比对:顺序实验/对照，已设对照为1
      tempDEG <- edgeR::topTags(lt, n = Inf,sort.by = 'logFC') 
      tempDEG <- as.data.frame(tempDEG)
      DEG_edgeR <- na.omit(tempDEG)
      DEG_edgeR$gene <- rownames(DEG_edgeR)
      
      select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.5,FDR < 0.05,gene %in% HEG)
      nrow(select_DEG_edgeR)
      if(nrow(select_DEG_edgeR) < 50){
        select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.2,PValue < 0.01)
      }
      
      if(nrow(select_DEG_edgeR) > 200){
        select_DEG_edgeR <- select_DEG_edgeR[1:200,]
      }
      
      markers_edgeR <- rownames(select_DEG_edgeR)
    }
    return(markers_edgeR)
    
  },mc.cores = threads)
  names(edgeR_marker_list) <- all_cts
  # sapply(edgeR_marker_list,length)
  
  # get neighbor_bg_genes for each ct
  neighbor_markers <- list()
  for(ct in all_cts){
    if(length(neighbor_celltypes[[ct]]) == 1){
      neighbor_markers[[ct]]$neighbor_markers <- global_markers[[ct]]$global_markers
      neighbor_markers[[ct]]$neighbor_bg_genes <- global_markers[[ct]]$global_bg_genes
    } else{
      neighbor_cts <- neighbor_celltypes[[ct]][-1]
      neighbor_bg_genes_tmp <- unlist(edgeR_marker_list[neighbor_cts],use.names = F) %>% unique()
      neighbor_markers[[ct]]$neighbor_markers <- edgeR_marker_list[[ct]]
      neighbor_markers[[ct]]$neighbor_bg_genes <- setdiff(neighbor_bg_genes_tmp,edgeR_marker_list[[ct]])
    }
  }
  
  return(neighbor_markers)
  
}







