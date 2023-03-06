########### GetSeedCandidates ###########
# get_seed_candidates


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


