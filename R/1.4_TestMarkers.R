########### TestMarkers ###########
# test_markers



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
  query_mtx <- NormalizeData(query_mtx,verbose = F)
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
  rm(scale.query_mtx)
  
  return(cell_meta)
}



