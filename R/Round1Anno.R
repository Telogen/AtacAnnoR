########### Round1Anno ###########
# get_cell_subtype
# get_all_markers

# get_feature_selected_pb_ref
# get_feature_selected_query_mtx
# get_cor_mtx
# get_kendall_pred
# get_cell_meta_with_true

# test_markers
# get_seed_barcode


#' Get cell subtypes from a pseudo bulk matrix
#'
#' Get cell subtypes from a pseudo bulk matrix
#'
#' @param pb.ref a pseudo bulk matrix whose rows are features and columns are cell types
#' @param min.cor the minimum cutoff of cell type correlation
#' @param verbose whether to display messages and plot, default is TRUE
#'
#' @return Returns a list of cell subtypes.
#' @export
#'
#' @importFrom pcaPP cor.fk
get_cell_subtype <- function(pb.ref,
                             min.cor = 0.6,
                             verbose = T) {
  COR <- pcaPP::cor.fk(pb.ref)
  cell_subtype_list <- apply(COR,1,function(row){
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
      ggplot2::theme_void() +
      ggplot2::labs(x = " ", y = " ") +
      ggplot2::theme(axis.text.y  = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,size = 10),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank()) +
      ggplot2::geom_text(ggplot2::aes(label = label),size = 3)
    print(p)
  }
  return(cell_subtype_list)
}



#' Get all markers
#'
#' Get global markers and neighbor specific markers from reference
#'
#' @param pb.ref a pseudo bulk matrix whose rows are features and columns are cell types
#' @param atlas which atlas to use, default is \code{'H'}
#'  \itemize{
#'   \item{H: }{use the Human Cell Atlas}
#'   \item{M: }{use the Mouse Cell Atlas}
#' }
#' @param cell.subtype the list of cell subtypes get from \code{get_cell_subtype()}
#' @param HEG.prop the proportion of highly expressed genes to first find, default is 0.2
#' @param marker.num the number of markers to find of each cell type
#'
#' @return Returns a list containing global markers and neighbor specific markers.
#' @export
#'
get_all_markers <- function(pb.ref,
                            atlas = "H",
                            cell.subtype,
                            HEG.prop = 0.2,
                            marker.num = 200) {

  # global markers
  if (atlas == "H") {
    utils::data(HCA, envir = environment())
    atlas <- HCA
  } else {
    utils::data(MCA, envir = environment())
    atlas <- MCA
  }
  pb_ref_tmp <- Seurat::NormalizeData(pb.ref, verbose = F) %>%
    BiocGenerics::t() %>%
    Seurat::ScaleData(do.center = F, verbose = F) %>%
    BiocGenerics::t()
  common_genes <- intersect(rownames(pb_ref_tmp), rownames(atlas))
  pb_all <- cbind(pb_ref_tmp[common_genes, ], atlas[common_genes, ])
  global_markers_df <- get_pb_markers(pb_all,
                                      ct.to.find.mk = colnames(pb_ref_tmp),
                                      each.ct.top = marker.num, HEG.prop = HEG.prop
  )
  global_marker_list <- lapply(colnames(pb_ref_tmp), function(ct) {
    global_markers <- global_markers_df[, ct]
    background_genes <- setdiff(rownames(pb_ref_tmp), global_markers) %>% sample(marker.num)
    markers <- list(
      global_markers = global_markers,
      background_genes = background_genes
    )
    return(markers)
  })
  names(global_marker_list) <- colnames(pb_ref_tmp)

  # neighbor markers
  neighbor_markers_list <- lapply(names(cell.subtype), function(ct) {
    each_cell_subtype_i <- cell.subtype[[`ct`]]
    if (length(each_cell_subtype_i) == 1) {
      pb_ref <- pb.ref
    } else {
      pb_ref <- pb.ref[, each_cell_subtype_i]
    }
    markers_df <- get_pb_markers(pb_ref, each.ct.top = marker.num, HEG.prop = HEG.prop)
    ct_markers <- markers_df[, ct]
    neighbor_markers <- markers_df[, setdiff(colnames(markers_df), ct)] %>%
      BiocGenerics::t() %>%
      as.vector() %>%
      unique()
    neighbor_markers <- neighbor_markers[1:marker.num]
    markers <- list(
      this_ct_markers = ct_markers,
      neighbor_markers = neighbor_markers
    )
    return(markers)
  })
  names(neighbor_markers_list) <- names(cell.subtype)

  all_markers <- list(
    global_markers = global_marker_list,
    neighbor_markers = neighbor_markers_list
  )
  return(all_markers)
}


#' Get feature selected pseudo bulk reference
#'
#' Get feature selected pseudo bulk reference matrix to calculate cell-cell type correlation
#'
#' @param pb.ref a pseudo bulk matrix whose rows are features and columns are cell types
#' @param HEG.prop the proportion of highly expressed genes to first find, default is 0.2,
#' if the matrix is real bulk matrix it is recommended to set as 0.8
#' @param each.ct.top the number of top marker genes to find for a cell type
#'
#' @return Returns a feature selected pseudo bulk reference matrix.
#' @export
#'
#' @examples
#' feature_selected_pb_ref <- get_feature_selected_pb_ref(pb_ref)
get_feature_selected_pb_ref <- function(pb.ref,
                                        HEG.prop = 0.2,
                                        each.ct.top = 200) {
  ref_selected_genes <- get_pb_markers(pb.ref, HEG.prop = HEG.prop, each.ct.top = each.ct.top) %>%
    unlist(use.names = F)
  pb_ref_selected <- pb.ref[ref_selected_genes, ]
  return(pb_ref_selected)
}


#' Get feature selected query matrix
#'
#' Get feature selected query matrix to calculate cell-cell type correlation
#'
#' @param query.mtx query gene activity matrix whose rows are genes and columns are cells
#' @param query.nmf.embedding a dimension reduced matrix whose rows are query cells and
#' columns are factors (meta-programs)
#'
#' @return Returns a feature selected query matrix.
#' @export
#'
#' @examples
#' get_feature_selected_query_mtx <- get_feature_selected_query_mtx(query_mtx, query_nmf_embedding)
get_feature_selected_query_mtx <- function(query.mtx,
                                           query.nmf.embedding) {
  clusters_tmp <- Seurat::FindNeighbors(query.nmf.embedding, verbose = F)$snn %>%
    Seurat::FindClusters(verbose = F)
  clusters_que <- as.character(clusters_tmp[, 1])
  names(clusters_que) <- rownames(clusters_tmp)
  pb_que <- get_pseudo_bulk_mtx(query.mtx, cell.type = clusters_que)
  que_selected_genes <- get_pb_markers(pb_que, HEG.prop = 0.2, each.ct.top = 200) %>%
    unlist(use.names = F)
  query_mtx_selected <- query.mtx[que_selected_genes, ]
  return(query_mtx_selected)
}


#' Get first round annotation correlation matrix
#'
#' Get first round annotation correlation matrix of each query cell to each reference cell type using
#' Kendall correlation coefficient
#'
#' @param feature.selected.pb.ref the feature selected pseudo bulk reference
#' @param feature.selected.query.mtx the feature selected query matrix
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a correlation matrix whose rows are query cells and columns are
#' reference cell types.
#' @export
#'
get_cor_mtx <- function(feature.selected.pb.ref,
                        feature.selected.query.mtx,
                        threads = 10,
                        verbose = T) {
  intersect_genes <- intersect(rownames(feature.selected.pb.ref), rownames(feature.selected.query.mtx))
  feature.selected.pb.ref <- feature.selected.pb.ref[intersect_genes, ]
  feature.selected.query.mtx <- feature.selected.query.mtx[intersect_genes, ]

  if (verbose == T) {
    message(paste0("Intersect genes number: ", length(intersect_genes)))
    message(paste0("Computing Kendall correlation coefficient using ", threads, " threads..."))
  }

  celltypes <- colnames(feature.selected.pb.ref)
  RESULT <- pbmcapply::pbmclapply(celltypes, function(celltype) {
    celltype_exp <- feature.selected.pb.ref[, celltype]
    celltype_cor <- apply(feature.selected.query.mtx, 2, pcaPP::cor.fk, y = celltype_exp)
    return(celltype_cor)
  }, mc.cores = threads)

  out <- as.data.frame(RESULT)
  colnames(out) <- celltypes
  out <- out[colnames(feature.selected.query.mtx), ]
  return(out)
}



#' Get first round annotation labels
#'
#' Get first round annotation labels from correlation matrix
#'
#' @param cor.mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a cell metadata whose rows are query cells and columns are first round labels.
#' @export
#'
#'
get_kendall_pred <- function(cor.mtx) {
  cell_meta <- data.frame(row.names = rownames(cor.mtx))
  cell_meta$kendall_pred <- colnames(cor.mtx)[apply(cor.mtx, 1, function(row) {
    ifelse(any(is.na(row)), c(1), which.max(row))
  })]
  return(cell_meta)
}


#' Add the true labels to cell metadata
#'
#' @param cell.meta a cell metadata
#' @param true.label a vector of true cell labels
#' @param cor.mtx the correlation matrix get from \code{get_cor_mtx()}
#'
#' @return Returns a new cell metadata with true labels.
#' @export
#'
#'
get_cell_meta_with_true <- function(cell.meta,
                                    true.label,
                                    cor.mtx = NULL) {
  if (!is.null(cor.mtx)) {
    get_rank <- function(row, names, ct.num, whose.rank) {
      names(row) <- names
      row_cor <- as.numeric(row[1:ct.num])
      names(row_cor) <- names[1:ct.num]
      sorted_row <- sort(row_cor, decreasing = T)
      rank <- which(names(sorted_row) == row[whose.rank])
      return(rank)
    }
    ct_num <- ncol(cor.mtx)
    cor.mtx$true <- true.label
    cell.meta$true_rank <- apply(cor.mtx, 1, get_rank, names = colnames(cor.mtx), ct_num, whose.rank = "true")
  }

  cell.meta$true <- true.label
  if ('kendall_pred' %in% colnames(cell.meta)){
    cell.meta$kendall_pred_booltrue <- (cell.meta$kendall_pred == true.label)
  }
  return(cell.meta)
}



#' Test markers of first round annotation
#'
#' Test markers of first round annotation
#'
#' @param query.mtx query gene activity matrix whose rows are genes and columns are cells
#' @param cell.meta a cell metadata
#' @param all.markers all markers get from \code{get_all_markers()}
#' @param threads the number of threads, default is 10
#' @param verbose whether to display messages, default is TRUE
#'
#' @return Returns a new cell metadata with columns 'GMSS','NMSS','NMFC'.
#' @export
#'
#'
test_markers <- function(query.mtx,
                         cell.meta,
                         which.label = 'kendall_pred',
                         all.markers,
                         threads = 10,
                         verbose = T) {
  used_genes <- unique(unlist(all.markers, use.names = F))
  scale.query.mtx <- Seurat::ScaleData(query.mtx[used_genes, ], do.center = F, verbose = F)

  if (verbose == T) {
    message(paste0("Start testing markers using ", threads, " threads..."))
  }

  barcodes <- row.names(cell.meta)
  RESULT <- pbmcapply::pbmclapply(barcodes, function(barcode) {
    # barcode <- barcodes[1]
    barcode_cellype <- as.character(cell.meta[barcode, which.label])
    # global markers test
    barcode.G.markers <- all.markers$global_markers[[barcode_cellype]]$global_markers
    barcode.G.markers.exp <- scale.query.mtx[barcode.G.markers, barcode]
    barcode.background.genes <- all.markers$global_markers[[barcode_cellype]]$background_genes
    barcode.background.genes.exp <- scale.query.mtx[barcode.background.genes, barcode]
    global.test <- stats::wilcox.test(x = barcode.G.markers.exp, y = barcode.background.genes.exp, alternative = "great")
    GMSS <- -log10(max(global.test$p.value, 1e-100))

    # neighbor markers test
    barcode.N.markers <- all.markers$neighbor_markers[[barcode_cellype]]$this_ct_markers
    barcode.N.markers.exp <- scale.query.mtx[barcode.N.markers, barcode]
    neighbor.N.markers <- all.markers$neighbor_markers[[barcode_cellype]]$neighbor_markers
    neighbor.N.markers.exp <- scale.query.mtx[neighbor.N.markers, barcode]
    neighbor.test <- stats::wilcox.test(x = barcode.N.markers.exp, y = neighbor.N.markers.exp, alternative = "great")
    NMSS <- -log10(max(neighbor.test$p.value, 1e-100))
    NMFC <- mean(barcode.N.markers.exp, na.rm = T) / mean(neighbor.N.markers.exp, na.rm = T)
    if (NMFC == "NaN") {
      NMFC <- 0
    }
    if (NMFC == Inf) {
      NMFC <- 5
    }

    out <- c(GMSS, NMSS, NMFC)
    return(out)
  }, mc.cores = threads)

  out <- matrix(unlist(RESULT), nrow = 3)
  cell.meta$GMSS <- out[1, ]
  cell.meta$NMSS <- out[2, ]
  cell.meta$NMFC <- out[3, ]
  rm(scale.query.mtx)

  return(cell.meta)
}




#' get_candidate_seed_barcodes
#'
#' @param cell_meta a cell metadata
#' @param threads number of threads
#'
#' @return the candidate seed barcodes
#' @export
#'
get_candidate_seed_barcodes <- function(cell_meta,threads = 10){

  cell_meta_li <- split(cell_meta,cell_meta$kendall_pred)
  names(cell_meta_li)
  seed_barcodes <- pbmcapply::pbmclapply(cell_meta_li,function(cell_meta_i){
    # cell_meta_i <- cell_meta_li$`Activated CD4 T cells`
    # GMSS
    hist(cell_meta_i$GMSS,100)
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
      mclustBIC <- mclust::mclustBIC
      if(nrow(cell_meta_i_mini) < 10000){
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
      cell_meta_i_mini$NMSS_classification <- as.character(mclust_out1$classification)

      # ggplot2::ggplot(cell_meta_i_mini, ggplot2::aes(x = NMSS, color = NMSS_classification)) +
      #   ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count), fill = NMSS_classification,alpha = 0.2)) +
      #   ggplot2::scale_y_continuous(expand = c(0,0)) +
      #   ggplot2::geom_vline(xintercept = final_NMSS_cutoff, linetype="dotted")

      # (ggplot(cell_meta_i_mini) +
      #     geom_point(aes(x = GMSS,y = NMSS,color = kendall_pred_booltrue)) +
      #     geom_vline(xintercept=2, linetype="dotted") +
      #     geom_hline(yintercept=final_NMSS_cutoff, linetype="dotted")) %>%
      #   aplot::insert_right((ggplot(cell_meta_i_mini,aes(x = NMSS,color = kendall_pred_booltrue)) +
      #                          geom_density(aes(y = after_stat(count),fill = kendall_pred_booltrue,alpha = 0.2)) +
      #                          scale_y_continuous(expand = c(0,0)) +
      #                          theme_void() +
      #                          theme(axis.title = element_blank(),
      #                                axis.text = element_blank(),
      #                                axis.ticks = element_blank()) +
      #                          coord_flip() +
      #                          NoLegend()),
      #                       0.2)
      # get_benchmark(cell_meta_i$true,cell_meta_i$kendall_pred)[1]
      # get_benchmark(cell_meta_i_mini$true,cell_meta_i_mini$kendall_pred)[1]
      # best_cluster_idx
      seed_cell_idx <- which(cell_meta_i_mini$NMSS >= final_NMSS_cutoff)
      # get_benchmark(cell_meta_i_mini[seed_cell_idx,]$true,
      #               cell_meta_i_mini[seed_cell_idx,]$kendall_pred)[1]
      length(seed_cell_idx)
      seed_cell_barcode <- rownames(cell_meta_i_mini)[seed_cell_idx]

    }
    return(seed_cell_barcode)
  },mc.cores = threads) %>% unlist(use.names = F)

  return(seed_barcodes)
}



