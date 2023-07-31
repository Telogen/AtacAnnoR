########### Main ###########
# RunAtacAnnoR
# RunAtacAnnoR_Signac
# RunAtacAnnoR_SnapATAC
# RunAtacAnnoR_ArchR



#' Run AtacAnnoR in a line
#'
#' AtacAnnoR main function.
#'
#' @param ref_mtx A reference gene expression counts matrix whose rows are genes and columns are cells/samples
#' @param ref_celltype A vector of cell types whose orders are same as the columns of `ref_mtx`
#' @param ref_type Which type of the reference matrix is
#'  \itemize{
#'   \item{'sc': }{The reference matrix is a single cell RNA-seq matrix}
#'   \item{'bulk': }{The reference matrix is a bulk RNA-seq matrix}
#' }
#' @param query_gene_activity A query gene activity matrix whose rows are genes and columns are cells
#' @param query_peak_counts A query scATAC-seq peak counts matrix whose rows are peaks and columns are cells
#' @param query_cell_embedding A query scATAC-seq cell embedding whose rows are cells and columns are factors, 
#' default is NULL, which means getting NMF embedding automatically
#' @param min_cor The minimum of correlation to define similar cell types
#' @param num_global_markers The number of global markers
#' @param num_neighbor_markers The number of neighbor markers
#' @param threads The number of threads
#' @param verbose Whether to display messages and plots
#' @param simple_output Whether to return a simple output
#'
#' @return If \code{simple_output} is set to TRUE, returns a vector of predicted labels.
#'   If \code{simple_output} is set to FALSE, returns a list object of all outputs.
#' @export
#'
RunAtacAnnoR <- function(ref_mtx, ref_celltype, ref_type = "sc",
                         query_gene_activity, query_peak_counts = NULL, query_cell_embedding = NULL,
                         min_cor = 0.6,num_global_markers = 200, num_neighbor_markers = 200,
                         threads = 10, verbose = TRUE, simple_output = TRUE) {
  if(ref_type == 'bulk'){
    get_global_markers <- get_global_markers_bulk
    get_neighbor_markers <- get_neighbor_markers_bulk
  } else{
    get_global_markers <- get_global_markers_sc
    get_neighbor_markers <- get_neighbor_markers_sc
  }
  
  # Pre-Processing
  if (verbose) {
    message("Pre-processing...")
  }
  pre_processing_mtxs <- pre_processing(ref_mtx = ref_mtx, query_mtx = query_gene_activity, verbose = verbose)
  ref_mtx <- pre_processing_mtxs$ref_mtx
  query_mtx <- pre_processing_mtxs$query_mtx
  rm(pre_processing_mtxs)
  gc()
  
  if (is.null(query_peak_counts) & is.null(query_cell_embedding)){
    stop("Please provide query_peak_counts or query_cell_embedding!")
  } else if (is.null(query_peak_counts)){
    if (verbose) {
      message("Using the given cell embedding")
    }
    query_nmf_embedding <- query_cell_embedding
  } else if (is.null(query_cell_embedding)){
    if (verbose) {
      message("Getting NMF embedding...")
    }
    query_nmf_embedding <- get_nmf_embedding(query_peak_counts, n_factors = 50, verbose = verbose)
  }
  
  ################### First Round Annotation ###################
  if (verbose) {
    message("First Round Annotation")
    message("Getting global and neighbor markers...")
  }
  global_markers <- get_global_markers(sc_counts_mtx = ref_mtx,labels = ref_celltype,max_marker = num_global_markers,
                                       threads = threads)
  sapply(global_markers,function(x){c(length(x[[1]]),length(x[[2]]))})
  neighbor_celltypes <- get_neighbor_celltypes(sc_count_mtx = ref_mtx,labels = ref_celltype,global_markers = global_markers,
                                               min_cor = min_cor,verbose = verbose)
  neighbor_markers <- get_neighbor_markers(sc_counts_mtx = ref_mtx,labels = ref_celltype,
                                           neighbor_celltypes = neighbor_celltypes,global_markers = global_markers,
                                           max_marker = num_neighbor_markers,threads = threads)
  
  if (verbose) {
    message("Getting candidate cell type labels...")
  }
  cor_mtx <- get_cor_mtx(sc_count_mtx = ref_mtx,labels = ref_celltype,query_mtx = query_mtx,
                         global_markers = global_markers,query_nmf_embedding = query_nmf_embedding,
                         threads = threads,verbose = verbose)
  cell_meta <- get_kendall_pred(cor_mtx)
  
  
  if (verbose) {
    message("Validating of candidate cell type labels...")
  }
  cell_meta <- test_markers(query_mtx,cell_meta,global_markers,neighbor_markers,
                            threads = threads,verbose = verbose)
  cell_meta <- get_seed_candidates(cell_meta, threads = threads)
  seed_candidate_meta <- cell_meta[which(cell_meta$is_seed_candidate == T),]
  
  ################### Second Round Annotation ###################
  if (verbose) {
    message("Cleaning seed cell candidates...")
  }
  cell_meta <- seed_cleaning(cell_meta,query_nmf_embedding)
  seed_meta <- cell_meta[which(cell_meta$is_seed == T),]
  
  if (verbose) {
    message("Final prediction using WKNN...")
  }
  cell_meta <- WKNN_predict(cell_meta,query_nmf_embedding,k = 10)
  
  
  # output
  if (simple_output) {
    out <- cell_meta$final_pred
  } else {
    out <- list(
      query_nmf_embedding = query_nmf_embedding,
      global_markers = global_markers,
      neighbor_markers = neighbor_markers,
      cor_mtx = cor_mtx,
      cell_meta = cell_meta
    )
    
  }
  return(out)
}





#' Run AtacAnnoR on a Seurat object in Signac pipeline
#'
#' Run AtacAnnoR on a Seurat object in Signac pipeline. 
#'
#' @param query_SeuratObj The query Seurat object
#' @param query_ga_assay The assay containing gene activity of the query Seurat object
#' @param query_peak_assay The assay containing peak counts of the query Seurat object
#' @param ref_SeuratObj The reference Seurat object
#' @param ref_assay The assay of the reference Seurat object to use
#' @param ref_ident The identity of the reference Seurat object which contains the reference cell type
#' @param ref_type Which type of the reference matrix is
#'  \itemize{
#'   \item{'sc': }{The reference matrix is a single cell RNS-seq matrix}
#'   \item{'bulk': }{The reference matrix is a bulk RNA-seq matrix}
#' }
#' @param min_cor The minimum of correlation to define similar cell types
#' @param num_global_markers The number of global markers
#' @param num_neighbor_markers The number of neighbor markers
#' @param threads The number of threads
#' @param verbose Whether to display messages and plots
#'
#' @return Returns a new query Seurat object with cell types predicted by AtacAnnoR restored in
#' \code{query$final_pred} and NMF embedding restored in \code{query[['nmf']]}.
#' @export
RunAtacAnnoR_Signac <- function(query_SeuratObj,query_ga_assay = 'ACTIVITY',query_peak_assay = 'ATAC',
                                ref_SeuratObj,ref_assay = 'RNA',ref_ident = 'celltype',ref_type = "sc",
                                min_cor = 0.6,num_global_markers = 200, num_neighbor_markers = 200,
                                threads = 10, verbose = TRUE){
  
  ref_mtx <- ref_SeuratObj[[`ref_assay`]]@counts
  ref_celltype <- ref_SeuratObj@meta.data[,ref_ident]
  
  query_gene_activity <- query_SeuratObj[[`query_ga_assay`]]@counts
  query_peak_counts <- query_SeuratObj[[`query_peak_assay`]]@counts
  
  pred <- RunAtacAnnoR(
    ref_mtx = ref_mtx, ref_celltype = ref_celltype, ref_type = ref_type,
    query_gene_activity = query_gene_activity, query_peak_counts = query_peak_counts, 
    min_cor = min_cor,num_global_markers = num_global_markers, num_neighbor_markers = num_neighbor_markers,
    threads = threads, verbose = verbose, simple_output = F)
  
  query_SeuratObj@meta.data <- cbind(query_SeuratObj@meta.data,pred$cell_meta)
  query_SeuratObj@reductions[['nmf']] <- Seurat::CreateDimReducObject(embeddings = pred$query_nmf_embedding,
                                                                      assay = query_peak_assay)
  return(query_SeuratObj)
}




#' Run AtacAnnoR on a snap object in SnapATAC pipeline
#'
#' Run AtacAnnoR on a snap object in SnapATAC pipeline. Note that before using this
#' function, the \code{gmat} and \code{pmat} should be stored in
#' \code{query_snapObj}.
#'
#' @param query_snapObj The query snap object
#' @param ref_mtx A reference gene expression matrix whose rows are genes and columns are cells/samples
#' @param ref_celltype A vector of cell types whose orders are same as the columns of `ref_mtx`
#' @param ref_type Which type of the reference matrix is
#'  \itemize{
#'   \item{'sc': }{The reference matrix is a single cell RNS-seq matrix}
#'   \item{'bulk': }{The reference matrix is a bulk RNA-seq matrix}
#' }
#' @param min_cor The minimum of correlation to define similar cell types
#' @param num_global_markers The number of global markers
#' @param num_neighbor_markers The number of neighbor markers
#' @param threads The number of threads
#' @param verbose Whether to display messages and plots
#'
#' @return Returns a new query snap object with cell types predicted by AtacAnnoR restored in
#' \code{query_snapObj@metaData$final_pred}.
#' @export
#'
RunAtacAnnoR_SnapATAC <- function(query_snapObj,ref_mtx, ref_celltype,ref_type = "sc",
                                  min_cor = 0.6,num_global_markers = 200, num_neighbor_markers = 200,
                                  threads = 10, verbose = TRUE){
  
  query_gene_activity <- BiocGenerics::t(query_snapObj@gmat)
  query_peak_counts <- BiocGenerics::t(query_snapObj@pmat)
  
  pred <- RunAtacAnnoR(
    ref_mtx = ref_mtx, ref_celltype = ref_celltype, ref_type = ref_type,
    query_gene_activity = query_gene_activity, query_peak_counts = query_peak_counts, 
    min_cor = min_cor,num_global_markers = num_global_markers, num_neighbor_markers = num_neighbor_markers,
    threads = threads, verbose = verbose, simple_output = F)
  
  query_snapObj@metaData <- cbind(query_snapObj@metaData,pred$cell_meta)
  return(query_snapObj)
}



#' Run AtacAnnoR on an ArchRProject in ArchR pipeline
#'
#' Run AtacAnnoR on an ArchRProject in ArchR pipeline Note that before using this
#' function, the \code{GeneScoreMatrix} and \code{PeakMatrix} should be stored in
#' \code{query_ArchRproj}, which can be checked by \code{getAvailableMatrices(query_ArchRproj)}.

#' @param query_ArchRproj The query ArchRProject
#' @param ref_mtx A reference gene expression matrix whose rows are genes and columns are cells/samples
#' @param ref_celltype A vector of cell types whose orders are same as the columns of `ref_mtx`
#' @param ref_type Which type of the reference matrix is
#'  \itemize{
#'   \item{'sc': }{The reference matrix is a single cell RNS-seq matrix}
#'   \item{'bulk': }{The reference matrix is a bulk RNA-seq matrix}
#' }
#' @param min_cor The minimum of correlation to define similar cell types
#' @param num_global_markers The number of global markers
#' @param num_neighbor_markers The number of neighbor markers
#' @param threads The number of threads
#' @param verbose Whether to display messages and plots
#'
#' @return Returns a new query ArchRProject with cell types predicted by AtacAnnoR restored in
#' \code{query_ArchRproj$final_pred}.
#' @export
#'
RunAtacAnnoR_ArchR <- function(query_ArchRproj,ref_mtx, ref_celltype,ref_type = "sc",
                               min_cor = 0.6,num_global_markers = 200, num_neighbor_markers = 200,
                               threads = 10, verbose = TRUE){
  
  AvailableMatrices <- ArchR::getAvailableMatrices(query_ArchRproj)
  if(min(c('GeneScoreMatrix','PeakMatrix') %in% AvailableMatrices) == 0){
    stop('GeneScoreMatrix and PeakMatrix must be present in query_ArchRproj!')
  }
  
  GeneScoreMatrix <- ArchR::getMatrixFromProject(query_ArchRproj,useMatrix = "GeneScoreMatrix",verbose = verbose)
  query_gene_activity <- GeneScoreMatrix@assays@data@listData$GeneScoreMatrix
  rownames(query_gene_activity) <- GeneScoreMatrix@elementMetadata@listData$name
  colnames(query_gene_activity) <- query_ArchRproj$cellNames

  PeakMatrix <- ArchR::getMatrixFromProject(query_ArchRproj,useMatrix = "PeakMatrix",verbose = verbose)
  query_peak_counts <- PeakMatrix@assays@data@listData$PeakMatrix
  rownames(query_peak_counts) <- PeakMatrix@elementMetadata@listData$name
  colnames(query_peak_counts) <- query_ArchRproj$cellNames
  
  pred <- RunAtacAnnoR(
    ref_mtx = ref_mtx, ref_celltype = ref_celltype, ref_type = ref_type,
    query_gene_activity = query_gene_activity, query_peak_counts = query_peak_counts, 
    min_cor = min_cor,num_global_markers = num_global_markers, num_neighbor_markers = num_neighbor_markers,
    threads = threads, verbose = verbose, simple_output = F)
  
  query_ArchRproj@cellColData <- cbind(query_ArchRproj@cellColData,pred$cell_meta)
  return(query_ArchRproj)
}



#' Run AtacAnnoR on a cell_data_set(cds) object in Cicero pipeline.
#'
#' Run AtacAnnoR on a cell_data_set(cds) object in Cicero pipeline. Note that before using this
#' function, gene activity should be calculate by Cicero (\code{https://cole-trapnell-lab.github.io/cicero-release/docs_m3/}).
#'
#' @param query_cds The query cds object
#' @param ref_mtx A reference gene expression matrix whose rows are genes and columns are cells/samples
#' @param ref_celltype A vector of cell types whose orders are same as the columns of `ref_mtx`
#' @param ref_type Which type of the reference matrix is
#'  \itemize{
#'   \item{'sc': }{The reference matrix is a single cell RNS-seq matrix}
#'   \item{'bulk': }{The reference matrix is a bulk RNA-seq matrix}
#' }
#' @param query_gene_activity A query gene activity matrix whose rows are genes and columns are cells
#' @param min_cor The minimum of correlation to define similar cell types
#' @param num_global_markers The number of global markers
#' @param num_neighbor_markers The number of neighbor markers
#' @param threads The number of threads
#' @param verbose Whether to display messages and plots
#'
#' @return Returns a new query cds object with cell types predicted by AtacAnnoR restored in
#' \code{query_cds$final_pred} and NMF embedding restored in \code{query_cds@int_colData$reducedDims$NMF}.
#' @export
RunAtacAnnoR_Cicero <- function(query_cds,ref_mtx, ref_celltype,ref_type = "sc",
                                min_cor = 0.6,num_global_markers = 200, num_neighbor_markers = 200,
                                query_gene_activity,threads = 10, verbose = TRUE){
  
  query_peak_counts <- query_cds@assays@data$counts
  
  pred <- RunAtacAnnoR(
    ref_mtx = ref_mtx, ref_celltype = ref_celltype, ref_type = ref_type,
    query_gene_activity = query_gene_activity, query_peak_counts = query_peak_counts, 
    min_cor = min_cor,num_global_markers = num_global_markers, num_neighbor_markers = num_neighbor_markers,
    threads = threads, verbose = verbose, simple_output = F)
  
  query_cds@colData <- cbind(query_cds@colData,pred$cell_meta)
  query_cds@int_colData$reducedDims$NMF <- pred$query_nmf_embedding
  return(query_cds)
}






