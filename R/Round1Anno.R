########### Round1Anno ###########



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








#' Get global markers for bulk RNA-seq gene counts matrix
#'
#' Get global markers for bulk RNA-seq gene counts matrix
#' 
#' @param sc_counts_mtx bulk RNA-seq counts matrix
#' @param labels bulk RNA-seq sample labels
#' @param max_marker maximum number of markers for each cell type
#' @param threads the number of threads, default is 10
#' @param return_raw whether return raw data, default is FALSE
#'
#' @return Returns a list of global markers for each cell type
#' @export
#'
get_global_markers_bulk <- function(sc_counts_mtx,labels,max_marker = 200,threads = 10,return_raw = F){
  
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





#' Get neighbor markers for bulk RNA-seq gene counts matrix
#'
#' Get neighbor markers for bulk RNA-seq gene counts matrix
#'
#' @param sc_counts_mtx bulk RNA-seq counts matrix
#' @param labels bulk RNA-seq sample labels
#' @param neighbor_celltypes neighbor celltypes got from `get_neighbor_celltypes()`
#' @param global_markers global markers got from `get_global_markers_bulk()`
#' @param max_marker maximum number of markers for each cell type
#' @param threads the number of threads, default is 10
#'
#' @return Returns a list of neighbor markers for each cell type
#' @export
#'
get_neighbor_markers_bulk <- function(sc_counts_mtx,labels,neighbor_celltypes, global_markers,max_marker = 200,threads = 10){
  
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
      keep.exprs <- edgeR::filterByExpr(dge) 
      dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
      dge <- edgeR::calcNormFactors(dge, method = 'TMM') 
      dge <- edgeR::estimateDisp(dge, design, robust=T) 
      fit <- edgeR::glmFit(dge, design, robust=T)
      lt <- edgeR::glmLRT(fit, contrast=c(-1,1))  
      tempDEG <- edgeR::topTags(lt, n = Inf,sort.by = 'logFC') 
      tempDEG <- as.data.frame(tempDEG)
      DEG_edgeR <- na.omit(tempDEG)
      DEG_edgeR$gene <- rownames(DEG_edgeR)
      
      select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.5,FDR < 0.05,gene %in% HEG)
      nrow(select_DEG_edgeR)
      if(nrow(select_DEG_edgeR) < 50){
        select_DEG_edgeR <- dplyr::filter(DEG_edgeR,logFC > 0.2,PValue < 0.01)
      }
      
      if(nrow(select_DEG_edgeR) > max_marker){
        select_DEG_edgeR <- select_DEG_edgeR[1:max_marker,]
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









#' Get correlation matrix
#'
#' Get correlation matrix of each query cell to each reference cell type using Kendall correlation coefficient
#'
#' @param sc_count_mtx reference gene counts matrix
#' @param labels reference cell type labels
#' @param query_mtx query_mtx query gene activity matrix
#' @param global_markers global markers got from `get_global_markers_sc()` or `get_global_markers_bulk()`
#' @param query_nmf_embedding query meta-program matrix
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a correlation matrix whose rows are query cells and columns are reference cell types
#' @export
#'
get_cor_mtx <- function(sc_count_mtx,labels,query_mtx,global_markers,query_nmf_embedding,threads = 10,verbose = T) {
  # ref
  if (verbose) {
    message("Processing reference...")
  }
  ref_selected_features <- unlist(global_markers) %>% unique()
  pb_ref_counts <- get_pseudo_bulk_mtx(sc_count_mtx,labels,mode = 'sum')
  
  # query
  if (verbose) {
    message("Processing query...")
  }
  query_clusters <- (Seurat::FindNeighbors(query_nmf_embedding, verbose = F)$snn %>%
                       Seurat::FindClusters(verbose = F))[,1] %>% as.character()
  query_pb <- get_pseudo_bulk_mtx(query_mtx, query_clusters,mode = 'sum')
  highly_expressed_genes <- apply(query_pb, 2, function(col) {
    sort(col, T)[1:ceiling(nrow(query_pb) * 0.2)] %>% names()
  }) %>% as.vector() %>% unique()
  query_pb.use <- query_pb[highly_expressed_genes, ]
  ideal_gene_mtx <- diag(ncol(query_pb.use))
  colnames(ideal_gene_mtx) <- colnames(query_pb.use)
  cos_mtx <- RcppML::cosine(t(query_pb.use), ideal_gene_mtx)
  dimnames(cos_mtx) <- list(rownames(query_pb.use),colnames(ideal_gene_mtx))
  query_selected_features <- apply(cos_mtx, 2, function(col) {
    sort(col, T)[1:200] %>% names()
  }) %>% as.data.frame() %>% unlist(use.names = F) %>% unique()
  # query_global_markers <- get_global_markers_sc(sc_counts_mtx = query_mtx,labels = query_clusters)
  # query_selected_features <- unlist(query_global_markers,use.names = F) %>% unique()
  
  
  # intersect_genes
  intersect_genes <- intersect(ref_selected_features, query_selected_features)
  feature_selected_pb_ref <- pb_ref_counts[intersect_genes, ]
  feature_selected_query_mtx <- query_mtx[intersect_genes, ]
  
  if (verbose) {
    message(paste0("Intersect genes number: ", length(intersect_genes)))
    message(paste0("Computing Kendall correlation coefficient using ", threads, " threads..."))
  }
  # correlation
  celltypes <- colnames(feature_selected_pb_ref)
  out <- pbmcapply::pbmclapply(celltypes, function(celltype) {
    celltype_exp <- feature_selected_pb_ref[, celltype]
    celltype_cor <- apply(feature_selected_query_mtx, 2, pcaPP::cor.fk, y = celltype_exp)
    return(celltype_cor)
  }, mc.cores = threads) %>% do.call(what = 'cbind')
  colnames(out) <- celltypes
  out <- as.data.frame(out)
  return(out)
}





#' Get first round annotation labels
#'
#' Get first round annotation labels from the correlation matrix
#'
#' @param cor_mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a cell metadata whose rows are query cells and columns are first round labels.
#' @export
#'
#'
get_kendall_pred <- function(cor_mtx) {
  cell_meta <- data.frame(row.names = rownames(cor_mtx))
  cell_meta$kendall_pred <- colnames(cor_mtx)[apply(cor_mtx, 1, function(row) {
    ifelse(any(is.na(row)), c(1), which.max(row))
  })]
  return(cell_meta)
}











#' Test markers
#'
#' Test whether the markers are of high activity in the predicted cells
#'
#' @param query_mtx query gene activity matrix whose rows are genes and columns are cells
#' @param cell_meta a cell metadata
#' @param global_markers global markers
#' @param neighbor_markers neighbor markers
#' @param which_label which label to test, default is 'kendall_pred'
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a new cell metadata with new columns 'GMSS','NMSS'
#' @export
#'
#'
test_markers <- function(query_mtx,cell_meta,global_markers,neighbor_markers,which_label = 'kendall_pred',threads = 10,verbose = T) {
  query_mtx <- Seurat::NormalizeData(query_mtx,verbose = F)
  used_genes <- unique(unlist(c(global_markers,neighbor_markers), use.names = F))
  scale.query_mtx <- Seurat::ScaleData(query_mtx[used_genes, ], do.center = F, verbose = F)
  # scale.query_mtx <- as.matrix(query_mtx[used_genes, ])
  
  if (verbose == T) {
    message(paste0("Start testing markers using ", threads, " threads..."))
  }
  
  barcodes <- rownames(cell_meta)
  RESULT <- pbmcapply::pbmclapply(barcodes, function(barcode) {
    # barcode <- barcodes[8]
    barcode_cellype <- as.character(cell_meta[barcode, which_label])
    # global markers test
    barcode.G.markers <- global_markers[[barcode_cellype]]$global_markers
    barcode.G.markers.exp <- scale.query_mtx[barcode.G.markers, barcode]
    barcode.GBG.genes <- global_markers[[barcode_cellype]]$global_bg_genes
    barcode.GBG.genes.exp <- scale.query_mtx[barcode.GBG.genes, barcode]
    global.test <- stats::wilcox.test(x = barcode.G.markers.exp, y = barcode.GBG.genes.exp, alternative = "great")
    GM.pvalue <- global.test$p.value
    -log10(GM.pvalue)
    
    # neighbor markers test
    barcode.N.markers <- neighbor_markers[[barcode_cellype]]$neighbor_markers
    barcode.N.markers.exp <- scale.query_mtx[barcode.N.markers, barcode]
    barcode.NBG.genes <- neighbor_markers[[barcode_cellype]]$neighbor_bg_genes
    barcode.NBG.genes.exp <- scale.query_mtx[barcode.NBG.genes, barcode]
    neighbor.test <- stats::wilcox.test(x = barcode.N.markers.exp, y = barcode.NBG.genes.exp, alternative = "great")
    NM.pvalue <- neighbor.test$p.value
    -log10(NM.pvalue)
    
    out <- c(GM.pvalue, NM.pvalue)
    return(out)
  }, mc.cores = threads)
  
  out <- matrix(unlist(RESULT), nrow = 2)
  cell_meta$GMSS <- stats::p.adjust(out[1, ],method = 'fdr') %>% sapply(function(x){-log10(max(x, 1e-100))})
  cell_meta$NMSS <- stats::p.adjust(out[2, ],method = 'fdr') %>% sapply(function(x){-log10(max(x, 1e-100))})
  
  return(cell_meta)
}








#' Get seed cell candidates
#' 
#' Determine GMSS cutoff and NMSSS cutoff for each cell type
#'
#' @param cell_meta a cell metadata
#' @param threads the number of threads, default is 10
#'
#' @return Returns a new cell_meta with column `GMSS_cutoff` and `NMSS_cutoff`
#' @export
#'
get_seed_candidates <- function(cell_meta,threads = 10){
  
  cell_meta_li <- split(cell_meta,cell_meta$kendall_pred)
  table(cell_meta$kendall_pred)
  names(cell_meta_li)
  new_cell_meta_li <- pbmcapply::pbmclapply(cell_meta_li,function(cell_meta_i){
    # cell_meta_i <- cell_meta_li[[3]]
    cell_meta_i$barcode <- rownames(cell_meta_i)
    cell_meta_i$GMSS_cutoff <- NULL
    cell_meta_i$NMSS_cutoff <- NULL
    # cell_meta_i$kendall_pred <- GML(cell_meta_i$kendall_pred)
    # cell_meta_i$true <- GML(cell_meta_i$true)
    # plot_cell_distribution(cell_meta_i)
    nrow(cell_meta_i)
    if(nrow(cell_meta_i) < 5){
      cell_meta_i$GMcluster <- '0'
      cell_meta_i$NMcluster <- '0'
      cell_meta_i$GMSS_cutoff <- Inf
      cell_meta_i$NMSS_cutoff <- Inf
    } else{
      # at least keep 20% cells after one filtering
      min_cell_num <- nrow(cell_meta_i)*0.2
      min_cell_num
      
      # GMSS
      suppressMessages(library(mclust))
      mclust_out1 <- mclust::Mclust(cell_meta_i$GMSS,G = 2:9,verbose = F)
      GMSS_cutoff_tmp <- sort(mclust_out1$parameters$mean,T)[1:2] %>% mean()
      GMSS_cutoff_candidates <- sort(c(sort(mclust_out1$parameters$mean,T)[-1],GMSS_cutoff_tmp),decreasing = T)
      for(GMSS_cutoff_i in GMSS_cutoff_candidates){
        kept_cell_num <- length(which(cell_meta_i$GMSS > GMSS_cutoff_i))
        if(kept_cell_num < min_cell_num){
          GMSS_cutoff <- GMSS_cutoff_i
        } else{
          GMSS_cutoff <- GMSS_cutoff_i
          break()
        }
      }
      GMSS_cutoff
      cell_meta_i$GMcluster <- as.character(mclust_out1$classification)
      # p1 <- ggplot(cell_meta_i, aes(x = GMSS, color = kendall_pred_booltrue)) +
      #   geom_density(aes(y = after_stat(count), fill = kendall_pred_booltrue,alpha = 0.2)) +
      #   scale_y_continuous(expand = c(0,0)) + 
      #   geom_vline(xintercept = mclust_out1$parameters$mean, linetype="dotted") + 
      #   geom_vline(xintercept = GMSS_cutoff, linetype=1) + 
      #   ggtitle(paste0(cell_meta_i$kendall_pred[1], ' kendall pred T/F distribution'))
      # p2 <- ggplot(cell_meta_i, aes(x = GMSS, color = GMcluster)) +
      #   geom_density(aes(y = after_stat(count), fill = GMcluster,alpha = 0.2)) +
      #   scale_y_continuous(expand = c(0,0)) + 
      #   geom_vline(xintercept=mclust_out1$parameters$mean, linetype="dotted") + 
      #   geom_vline(xintercept = GMSS_cutoff, linetype=1) + 
      #   ggtitle(paste0(cell_meta_i$kendall_pred[1], ' GMM cluster distribution'))
      # aplot::insert_bottom(p1,p2)
      cell_meta_i_GMSS <- dplyr::filter(cell_meta_i,GMSS > GMSS_cutoff)
      # acc1 <- get_benchmark(GML(cell_meta_i_GMSS$true),GML(cell_meta_i_GMSS$kendall_pred))[1]
      # print(paste0('GMSS cluster number: ',length(mclust_out1$parameters$mean)))
      # print(paste0('after GMSS cells: ',nrow(cell_meta_i_GMSS),'; after GMSS accuracy: ',round(acc1,4)))
      
      # NMSS
      mclust_out2 <- mclust::Mclust(cell_meta_i$NMSS,G = 2:9,verbose = F)
      NMSS_cutoff_tmp <- sort(mclust_out2$parameters$mean,T)[1:2] %>% mean()
      NMSS_cutoff_candidates <- sort(c(sort(mclust_out2$parameters$mean,T)[-1],NMSS_cutoff_tmp),decreasing = T)
      for(NMSS_cutoff_i in NMSS_cutoff_candidates){
        kept_cell_num <- length(which(cell_meta_i$NMSS > NMSS_cutoff_i))
        if(kept_cell_num < min_cell_num){
          NMSS_cutoff <- NMSS_cutoff_i
        } else{
          NMSS_cutoff <- NMSS_cutoff_i
          break()
        }
      }
      NMSS_cutoff
      cell_meta_i$NMcluster <- as.character(mclust_out2$classification)
      # p3 <- ggplot(cell_meta_i, aes(x = NMSS, color = kendall_pred_booltrue)) +
      #   geom_density(aes(y = after_stat(count), fill = kendall_pred_booltrue,alpha = 0.2)) +
      #   scale_y_continuous(expand = c(0,0)) + 
      #   geom_vline(xintercept = mclust_out2$parameters$mean, linetype="dotted") + 
      #   geom_vline(xintercept = NMSS_cutoff, linetype=1) + 
      #   ggtitle(paste0(cell_meta_i$kendall_pred[1], ' kendall pred T/F distribution'))
      # p4 <- ggplot(cell_meta_i, aes(x = NMSS, color = NMcluster)) +
      #   geom_density(aes(y = after_stat(count), fill = NMcluster,alpha = 0.2)) +
      #   scale_y_continuous(expand = c(0,0)) + 
      #   geom_vline(xintercept=mclust_out2$parameters$mean, linetype="dotted") + 
      #   geom_vline(xintercept = NMSS_cutoff, linetype=1) + 
      #   ggtitle(paste0(cell_meta_i$kendall_pred[1], ' GMM cluster distribution'))
      # aplot::insert_bottom(p3,p4)
      cell_meta_i_NMSS <- dplyr::filter(cell_meta_i_GMSS,NMSS > NMSS_cutoff)
      # acc2 <- get_benchmark(GML(cell_meta_i_NMSS$true),GML(cell_meta_i_NMSS$kendall_pred))[1]
      # print(paste0('NMSS cluster number: ',length(mclust_out2$parameters$mean)))
      # print(paste0('after NMSS cells: ',nrow(cell_meta_i_NMSS),'; after NMSS accuracy: ',round(acc2,4)))
      
      cell_meta_i$GMSS_cutoff <- GMSS_cutoff
      cell_meta_i$NMSS_cutoff <- NMSS_cutoff
    }
    
    return(cell_meta_i)
    # seed_cell_barcode <- rownames(cell_meta_i_NMSS)
    # return(seed_cell_barcode)
  },mc.cores = threads)
  
  new_cell_meta <- do.call(what = 'rbind',new_cell_meta_li)
  rownames(new_cell_meta) <- new_cell_meta$barcode
  cell_meta <- new_cell_meta[rownames(cell_meta),]
  cell_meta$barcode <- NULL
  cell_meta$is_seed_candidate <- FALSE
  cell_meta[which(cell_meta$GMSS > cell_meta$GMSS_cutoff & cell_meta$NMSS > cell_meta$NMSS_cutoff),]$is_seed_candidate <- TRUE
  
  # seed_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
  # get_benchmark(GML(cell_meta$true),GML(cell_meta$kendall_pred))
  # get_benchmark(GML(seed_meta$true),GML(seed_meta$kendall_pred))
  # nrow(seed_meta)/nrow(cell_meta)
  # get_each_recall(seed_meta)
  
  return(cell_meta)
}






