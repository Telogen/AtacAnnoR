########### Round1Anno ###########
# get_pseudo_bulk_mtx
# get_global_markers_pb
# get_neighbor_celltypes
# get_neighbor_markers_pb

# get_cor_mtx
# get_kendall_pred
# get_cell_meta_with_true

# test_markers
# get_seed_barcode



#' Get the pseudo bulk matrix
#'
#' Get the pseudo bulk matrix from a single cell matrix
#'
#' @param sc_mtx a single cell matrix whose rows are features and columns are cells
#' @param labels a vector of cell types whose orders are same as the columns of sc_mtx
#' @param mode the method to aggragate cells, default is \code{'mean'}
#'  \itemize{
#'   \item{'mean': }{get the mean of each cell type}
#'   \item{'sum': }{get the sum of each cell type}
#'   \item{'median': }{get the median of each cell type}
#' }
#'
#' @return Returns a pseudo bulk matrix whose rows are features and columns are cell types.
#' @export
#'
#' @examples
#' pb_ref <- get_pseudo_bulk_mtx(sc_mtx = ref_mtx, labels = REF_CELLTYPE)
get_pseudo_bulk_mtx <- function(sc_mtx, labels, mode = "mean") {
  Seurat_obj = SeuratObject::CreateSeuratObject(sc_mtx,meta.data = data.frame(group = labels,row.names = colnames(sc_mtx)))
  if (mode == "mean") {
    out_mtx <- Seurat::AverageExpression(Seurat_obj,group.by = 'group',slot = 'counts')[[1]]
  } else{
    out_mtx <- Seurat::AggregateExpression(Seurat_obj,group.by = 'group',slot = 'counts')[[1]]
  }
  return(out_mtx)
}




#' get_metacell_exp
#'
#' @param count_mtx todo
#' @param labels todo
#' @param meta_cell_num todo
#' @param mode todo
#' @param return_Seurat todo
#'
#' @return todo
#' @export
#'
#' @examples todo
get_metacell_exp <- function(count_mtx,labels,meta_cell_num = 10,mode = 'sum',return_Seurat = T){
  barcodes_list <- split(colnames(count_mtx),labels)
  meta_cell_exp <- lapply(names(barcodes_list),function(ct){
    ct_barcodes <- barcodes_list[[ct]]
    count_mtx_ct <- count_mtx[,ct_barcodes]
    meta_cell_labels <- rep(1:meta_cell_num,length.out = length(ct_barcodes)) %>%
      paste0(ct,'=',.) %>% sample()
    ct_meta_cell_exp <- get_pseudo_bulk_mtx(sc_mtx = count_mtx_ct,labels = meta_cell_labels,mode = mode)
    return(ct_meta_cell_exp)
  }) %>% do.call(what = 'cbind')
  if(return_Seurat){
    Seurat_obj <- CreateSeuratObject(meta_cell_exp,names.delim = '=')
    Seurat_obj <- NormalizeData(Seurat_obj,verbose = F)
    return(Seurat_obj)
  } else{
    return(meta_cell_exp)
  }
}




#' get_global_markers_pb
#'
#' @param metacell_Seurat todo
#' @param threads todo
#'
#' @return todo
#' @export
#'
#' @examples todo
get_global_markers_pb <- function(metacell_Seurat,threads = 10){
  pb_counts <- get_pseudo_bulk_mtx(metacell_Seurat@assays$RNA@counts,metacell_Seurat$orig.ident,mode = 'sum')
  all_cts <- levels(metacell_Seurat$orig.ident)
  Idents(metacell_Seurat) <- 'orig.ident'

  global_markers <- pbmcapply::pbmclapply(all_cts,function(ct){
    # ct <- 'CD4_TE'
    ct_HEG <- pb_counts[,ct] %>% sort(T) %>% names() %>% head(.,ceiling(length(.)*0.8))
    ct_global_mak_df <- FindMarkers(metacell_Seurat,features = ct_HEG,ident.1 = ct,only.pos = T,verbose = F,
                                    logfc.threshold = 0.5,min.cells.group = 1)
    ct_global_mak_df$fdr <- stats::p.adjust(ct_global_mak_df$p_val,method = 'fdr')
    if(length(which(ct_global_mak_df$fdr < 0.01)) > 10){
      ct_global_mak_df <- dplyr::filter(ct_global_mak_df,fdr < 0.01)
    }
    ct_global_markers <- rownames(ct_global_mak_df)
    ct_global_bg_genes <- sample(rownames(pb_counts),1000)

    ct_markers <- list(global_markers = ct_global_markers,
                       global_bg_genes = ct_global_bg_genes)
    return(ct_markers)
  },mc.cores = threads)
  names(global_markers) <- all_cts

  return(global_markers)
}



#' Get cell neighbors from a pseudo bulk matrix
#'
#' Get cell neighbors from a pseudo bulk matrix
#'
#' @param pb.ref a pseudo bulk matrix whose rows are features and columns are cell types
#' @param min.cor the minimum cutoff of cell type correlation
#' @param verbose whether to display messages and plot, default is TRUE
#'
#' @return Returns a list of cell subtypes.
#' @export
#'
#' @importFrom pcaPP cor.fk
get_neighbor_celltypes <- function(metacell_Seurat,global_markers,min.cor = 0.6,verbose = T) {
  pb_counts <- get_pseudo_bulk_mtx(metacell_Seurat@assays$RNA@counts,metacell_Seurat$orig.ident,mode = 'sum')
  all_global_markers <- unlist(lapply(global_markers,function(x){x$global_markers}),use.names = F) %>% unique()
  pb_data_selected <- pb_counts[all_global_markers,]
  COR <- pcaPP::cor.fk(pb_data_selected)
  neighbor_celltypes_list <- apply(COR,1,function(row){
    names(which(sort(row, decreasing = T) > min.cor))
  })
  if (verbose == T) {
    ord <- hclust(dist(COR, method = "euclidean"), method = "ward.D" )$order
    data <- as.data.frame(as.table(COR[ord,ord]))
    data$label <- round((data$Freq*100),0)
    data$label[which(data$label < (min.cor*100))] <- ' '
    p <- ggplot2::ggplot(data, ggplot2::aes(Var1, Var2, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(name = "cor",low = "blue", mid = 'yellow', high = "red")  +
      ggplot2::theme_bw() +
      ggplot2::labs(x = " ", y = " ") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line =  ggplot2::element_blank()) +
      ggplot2::geom_text( ggplot2::aes(label = label),size = 3)
    print(p)
  }
  return(neighbor_celltypes_list)
}



#' get_neighbor_markers_pb
#'
#' @param metacell_Seurat todo
#' @param neighbor_celltypes todo
#' @param global_markers todo
#' @param threads todo
#'
#' @return todo
#' @export
#'
#' @examples todo
get_neighbor_markers_pb <- function(metacell_Seurat,neighbor_celltypes,global_markers,threads = 10){
  pb_counts <- get_pseudo_bulk_mtx(metacell_Seurat@assays$RNA@counts,metacell_Seurat$orig.ident,mode = 'sum')
  all_cts <- levels(metacell_Seurat$orig.ident)
  Idents(metacell_Seurat) <- 'orig.ident'
  all_global_markers <- unlist(lapply(global_markers,function(x){x$global_markers}),use.names = F) %>% unique()

  neighbor_markers <- pbmcapply::pbmclapply(all_cts,function(ct){
    # ct <- 'CD4_TE'
    ct_HEG <- pb_counts[,ct] %>% sort(T) %>% names() %>% head(.,ceiling(length(.)*0.8))
    if(length(neighbor_celltypes[[ct]]) == 1){
      ct_neighbor_markers <- global_markers[[ct]]$global_markers
      ct_neighbor_bg_genes <- setdiff(all_global_markers,ct_neighbor_markers)
    } else{
      metacell_Seurat_sub <- subset(metacell_Seurat,cells = which(metacell_Seurat$orig.ident %in% neighbor_celltypes[[ct]]))
      ct_all_neighbor_mak_df <- FindAllMarkers(metacell_Seurat_sub,only.pos = T,verbose = F,
                                               logfc.threshold = 0.5,min.cells.group = 1,return.thresh = 0.1)
      ct_all_neighbor_mak_df$fdr <- stats::p.adjust(ct_all_neighbor_mak_df$p_val,method = 'fdr')
      if(min(table(dplyr::filter(ct_all_neighbor_mak_df,fdr < 0.01)$cluster)) > 10){
        ct_all_neighbor_mak_df <- dplyr::filter(ct_all_neighbor_mak_df,fdr < 0.01)
      }
      ct_neighbor_markers_tmp <- dplyr::filter(ct_all_neighbor_mak_df,cluster == ct)$gene
      ct_neighbor_bg_genes_tmp <- dplyr::filter(ct_all_neighbor_mak_df,cluster != ct)$gene %>% unique()
      ct_neighbor_markers <- setdiff(ct_neighbor_markers_tmp,ct_neighbor_bg_genes_tmp)
      ct_neighbor_bg_genes <- setdiff(ct_neighbor_bg_genes_tmp,ct_neighbor_markers_tmp)
    }

    ct_markers <- list(neighbor_markers = ct_neighbor_markers,
                       neighbor_bg_genes = ct_neighbor_bg_genes)
    return(ct_markers)
  },mc.cores = threads)
  names(neighbor_markers) <- all_cts
  return(neighbor_markers)
}






#' Get feature selected pseudo bulk reference
#'
#' Get first round annotation correlation matrix
#'
#' Get first round annotation correlation matrix of each query cell to each reference cell type using
#' Kendall correlation coefficient
#'
#' @param ref_metacell_Seurat ref_metacell_Seurat
#' @param query_mtx query_mtx
#' @param global_markers global_markers
#' @param query_mtx query_mtx
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a correlation matrix whose rows are query cells and columns are
#' reference cell types.
#' @export
#'
#' @examples
#' cor_mtx <- get_cor_mtx(feature_selected_pb_ref, feature_selected_query_mtx)
get_cor_mtx <- function(ref_metacell_Seurat,query_mtx,global_markers,query_nmf_embedding,threads = 10,verbose = T) {
  # ref
  if (verbose == T) {
    message("Processing reference...")
  }
  all_ref_global_markers <- unlist(lapply(global_markers,function(x){x$global_markers}),use.names = F) %>% unique()
  pb_ref_counts <- get_pseudo_bulk_mtx(ref_metacell_Seurat@assays$RNA@counts,ref_metacell_Seurat$orig.ident,mode = 'sum')
  # query
  if (verbose == T) {
    message("Processing query...")
  }
  query_clusters <- (Seurat::FindNeighbors(query_nmf_embedding, verbose = F)$snn %>%
                       Seurat::FindClusters(verbose = F))[,1]
  query_metacell_Seurat <- get_metacell_exp(query_mtx,query_clusters)
  query_global_markers <- get_global_markers_pb(query_metacell_Seurat,threads = threads)
  all_query_global_markers <- unlist(lapply(query_global_markers,function(x){x$global_markers}),use.names = F) %>% unique()
  # intersect_genes
  intersect_genes <- intersect(all_ref_global_markers, all_query_global_markers)
  feature.selected.pb.ref <- pb_ref_counts[intersect_genes, ]
  feature.selected.query_mtx <- query_mtx[intersect_genes, ]

  if (verbose == T) {
    message(paste0("Intersect genes number: ", length(intersect_genes)))
    message(paste0("Computing Kendall correlation coefficient using ", threads, " threads..."))
  }
  # correlation
  celltypes <- colnames(feature.selected.pb.ref)
  out <- pbmcapply::pbmclapply(celltypes, function(celltype) {
    celltype_exp <- feature.selected.pb.ref[, celltype]
    celltype_cor <- apply(feature.selected.query_mtx, 2, pcaPP::cor.fk, y = celltype_exp)
    return(celltype_cor)
  }, mc.cores = threads) %>% do.call(what = 'cbind')
  colnames(out) <- celltypes
  out <- as.data.frame(out)
  return(out)
}





#' Get first round annotation labels
#'
#' Get first round annotation labels from correlation matrix
#'
#' @param cor_mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a cell metadata whose rows are query cells and columns are first round labels.
#' @export
#'
#' @examples
#' cell_meta <- get_kendall_pred(cor_mtx)
#'
get_kendall_pred <- function(cor_mtx) {
  cell_meta <- data.frame(row.names = rownames(cor_mtx))
  cell_meta$kendall_pred <- colnames(cor_mtx)[apply(cor_mtx, 1, function(row) {
    ifelse(any(is.na(row)), c(1), which.max(row))
  })]
  return(cell_meta)
}


#' Add the true labels to cell metadata
#'
#' @param cell_meta a cell metadata
#' @param true_labels a vector of true cell labels
#' @param cor_mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a new cell metadata with true labels.
#' @export
#'
#' @examples
#' cell_meta <- get_cell_meta_with_true(cell_meta, true_labels = TRUE_LABEL, cor_mtx)
#'
get_cell_meta_with_true <- function(cell_meta,true_labels,cor_mtx = NULL) {
  if (!is.null(cor_mtx)) {
    get_rank <- function(row, names, ct.num, whose.rank) {
      names(row) <- names
      row_cor <- as.numeric(row[1:ct.num])
      names(row_cor) <- names[1:ct.num]
      sorted_row <- sort(row_cor, decreasing = T)
      rank <- which(names(sorted_row) == row[whose.rank])
      return(rank)
    }
    ct_num <- ncol(cor_mtx)
    cor_mtx$true <- true_labels
    cell_meta$true_rank <- apply(cor_mtx, 1, get_rank, names = colnames(cor_mtx), ct_num, whose.rank = "true")
  }

  cell_meta$true <- true_labels
  if ('kendall_pred' %in% colnames(cell_meta)){
    cell_meta$kendall_pred_booltrue <- (cell_meta$kendall_pred == true_labels)
  }
  return(cell_meta)
}



#' Test markers of first round annotation
#'
#' Test markers of first round annotation
#'
#' @param query_mtx query gene activity matrix whose rows are genes and columns are cells
#' @param cell_meta a cell metadata
#' @param global_markers global_markers
#' @param neighbor_markers neighbor_markers
#' @param which_label which_label
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a new cell metadata with columns 'GMSS','NMSS','NMFC'.
#' @export
#'
#'
test_markers <- function(query_mtx,cell_meta,global_markers,neighbor_markers,which_label = 'kendall_pred',threads = 10,verbose = T) {
  used_genes <- unique(unlist(c(global_markers,neighbor_markers), use.names = F))
  scale.query_mtx <- Seurat::ScaleData(query_mtx[used_genes, ], do.center = F, verbose = F)
  # scale.query_mtx <- as.matrix(query_mtx[used_genes, ])

  if (verbose == T) {
    message(paste0("Start testing markers using ", threads, " threads..."))
  }

  barcodes <- row.names(cell_meta)
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

    # neighbor markers test
    barcode.N.markers <- neighbor_markers[[barcode_cellype]]$neighbor_markers
    barcode.N.markers.exp <- scale.query_mtx[barcode.N.markers, barcode]
    barcode.NBG.genes <- neighbor_markers[[barcode_cellype]]$neighbor_bg_genes
    barcode.NBG.genes.exp <- scale.query_mtx[barcode.NBG.genes, barcode]
    neighbor.test <- stats::wilcox.test(x = barcode.N.markers.exp, y = barcode.NBG.genes.exp, alternative = "great")
    NM.pvalue <- neighbor.test$p.value

    out <- c(GM.pvalue, NM.pvalue)
    return(out)
  }, mc.cores = threads)

  out <- matrix(unlist(RESULT), nrow = 2)
  cell_meta$GMSS <- stats::p.adjust(out[1, ],method = 'fdr') %>% sapply(function(x){-log10(max(x, 1e-100))})
  cell_meta$NMSS <- stats::p.adjust(out[2, ],method = 'fdr') %>% sapply(function(x){-log10(max(x, 1e-100))})
  rm(scale.query_mtx)

  return(cell_meta)
}





#' get_seed_barcode
#'
#' @param cell_meta
#'
#' @return the seed barcodes
#' @export
#'
get_seed_candidates <- function(cell_meta){

  cell_meta_li <- split(cell_meta,cell_meta$kendall_pred)
  names(cell_meta_li)
  suppressMessages(library(mclust))
  seed_barcodes <- pbmcapply::pbmclapply(cell_meta_li,function(cell_meta_i){
    # cell_meta_i <- cell_meta_li[[8]]
    # GMSS
    plot_cell_distribution(cell_meta_i)
    hist(cell_meta_i$GMSS,100)
    plot_cell_distribution(cell_meta_i)
    get_quantile <- function(data,x){
      which.min(abs(sort(data) - x))/length(data)
    }
    get_quantile(cell_meta_i$GMSS,2)
    # ggplot(cell_meta_i, aes(x = GMSS, color = kendall_pred_booltrue)) +
    #   geom_density(aes(y = after_stat(count), fill = kendall_pred_booltrue,alpha = 0.2)) +
    #   scale_y_continuous(expand = c(0,0))
    nrow(cell_meta_i)
    if (nrow(cell_meta_i) < 200){
      cell_meta_i_mini <- cell_meta_i[which(cell_meta_i$GMSS > 1),]
    } else{
      cell_meta_i_mini <- cell_meta_i[which(cell_meta_i$GMSS > 2),]
    }

    # NMSS,NMFC
    nrow(cell_meta_i_mini)
    if(nrow(cell_meta_i_mini) < 3){
      seed_cell_barcode <- NULL
    } else{
      if(nrow(cell_meta_i_mini) < 100000){
        mclust_out1 <- mclust::Mclust(cell_meta_i_mini$NMSS,G = 2,verbose = F)
      } else{
        mclust_out1 <- mclust::Mclust(cell_meta_i_mini$NMSS,G = 10,verbose = F)
      }
      NMSS_cutoffs <- sort(mclust_out1$parameters$mean,decreasing = T)
      NMSS_cutoffs
      if(nrow(cell_meta_i_mini) < 200){
        min.number <- nrow(cell_meta_i_mini)/2
      } else{
        min.number <- 100
      }
      min.number
      for(NMSS_cutoff in NMSS_cutoffs){
        if(length(which(cell_meta_i_mini$NMSS >= NMSS_cutoff)) < min.number){
          final_NMSS_cutoff <- NMSS_cutoff
        } else{
          final_NMSS_cutoff <- NMSS_cutoff
          break()
        }
      }
      final_NMSS_cutoff
      seed_cell_idx <- which(cell_meta_i_mini$NMSS >= final_NMSS_cutoff)
      get_benchmark(cell_meta_i_mini[seed_cell_idx,]$true,
                    cell_meta_i_mini[seed_cell_idx,]$kendall_pred)[1]
      length(seed_cell_idx)
      seed_cell_barcode <- rownames(cell_meta_i_mini)[seed_cell_idx]
    }
    return(seed_cell_barcode)
  },mc.cores = 10) %>% unlist(use.names = F)

  return(seed_barcodes)
}







