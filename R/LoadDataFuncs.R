# load_data_funcs

GET_SeuratObj_RNA <- function(dataset){
  if(dataset == '10X-Multiome'){
    SeuratObj_RNA <- SeuratObj_RNA.1
  } else if(dataset == 'SNARE-seq'){
    SeuratObj_RNA <- SeuratObj_RNA.2
  } else if(dataset == 'SHARE-seq'){
    SeuratObj_RNA <- SeuratObj_RNA.3
  } else if(dataset == 'BloodDev'){
    SeuratObj_RNA <- SeuratObj_RNA.4
  } else if(dataset == 'Kidney'){
    SeuratObj_RNA <- SeuratObj_RNA.5
  } else if(dataset == 'Brain'){
    SeuratObj_RNA <- SeuratObj_RNA.6
  } else if(dataset == 'Zhu'){
    SeuratObj_RNA <- SeuratObj_RNA.7
  } else if(dataset == 'Wilk'){
    SeuratObj_RNA <- SeuratObj_RNA.8
  } else if(dataset == 'Stephenson'){
    SeuratObj_RNA <- SeuratObj_RNA.9
  } else if(dataset == 'Hao'){
    SeuratObj_RNA <- SeuratObj_RNA.10
  } else if(dataset == 'Monaco'){
    SeuratObj_RNA <- SeuratObj_RNA.11
  }
  return(SeuratObj_RNA)
}


GET_SeuratObj_ATAC <- function(dataset){
  if(dataset %in% c('10X-Multiome','Zhu','Wilk','Stephenson','Hao','Monaco')){
    SeuratObj_ATAC <- SeuratObj_ATAC.1
  } else if(dataset == 'SNARE-seq'){
    SeuratObj_ATAC <- SeuratObj_ATAC.2
  } else if(dataset == 'SHARE-seq'){
    SeuratObj_ATAC <- SeuratObj_ATAC.3
  } else if(dataset == 'BloodDev'){
    SeuratObj_ATAC <- SeuratObj_ATAC.4
  } else if(dataset == 'Kidney'){
    SeuratObj_ATAC <- SeuratObj_ATAC.5
  } else if(dataset == 'Brain'){
    SeuratObj_ATAC <- SeuratObj_ATAC.6
  }
  return(SeuratObj_ATAC)
}

GET_GML <- function(dataset){
  if(dataset %in% c('10X-Multiome','SNARE-seq','BloodDev','Brain')){
    GML <- function(x){return(x)}
  } else if(dataset == 'SHARE-seq'){
    GML <- get_merged_labels_dataset_SHAREseq
  } else if(dataset == 'Kidney'){
    GML <- get_merged_labels_Kidney
  } else if(dataset == 'Zhu'){
    GML <- get_merged_labels_Zhu
  } else if(dataset == 'Wilk'){
    GML <- get_merged_labels_Wilk
  } else if(dataset == 'Stephenson'){
    GML <- get_merged_labels_Stephenson
  } else if(dataset == 'Hao'){
    GML <- get_merged_labels_Hao
  } else if(dataset == 'Monaco'){
    GML <- get_merged_labels_Monaco
  }
  return(GML)
}




