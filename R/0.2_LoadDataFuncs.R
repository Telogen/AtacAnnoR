# load_data_funcs

GET_SeuratObj_RNA <- function(dataset,ls){
  if(dataset == '10X-Multiome'){
    if('SeuratObj_RNA.1' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.1
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/10X-Multiome/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'SNARE-seq'){
    if('SeuratObj_RNA.2' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.2
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/SNARE-seq/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'SHARE-seq'){
    if('SeuratObj_RNA.3' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.3
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/SHARE-seq/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'BloodDev'){
    if('SeuratObj_RNA.4' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.4
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/BloodDev/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Kidney'){
    if('SeuratObj_RNA.5' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.5
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Kidney/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Brain'){
    if('SeuratObj_RNA.6' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.6
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Brain/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Zhu'){
    if('SeuratObj_RNA.7' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.7
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Zhu/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Wilk'){
    if('SeuratObj_RNA.8' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.8
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Wilk/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Stephenson'){
    if('SeuratObj_RNA.9' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.9
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Stephenson/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Hao'){
    if('SeuratObj_RNA.10' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.10
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Hao/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'Monaco'){
    if('SeuratObj_RNA.11' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.11
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Monaco/SeuratObj/SeuratObj_RNA.RDS')
    }
  } else if(dataset == 'ISSAAC-seq'){
    if('SeuratObj_RNA.12' %in% ls){
      SeuratObj_RNA <- SeuratObj_RNA.12
    } else{
      SeuratObj_RNA <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/ISSAAC-seq/SeuratObj/SeuratObj_RNA.RDS')
    }
  }
  return(SeuratObj_RNA)
}


GET_SeuratObj_ATAC <- function(dataset,ls){
  if(dataset %in% c('10X-Multiome','Zhu','Wilk','Stephenson','Hao','Monaco')){
    if('SeuratObj_ATAC.1' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.1
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/10X-Multiome/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'SNARE-seq'){
    if('SeuratObj_ATAC.2' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.2
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/SNARE-seq/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'SHARE-seq'){
    if('SeuratObj_ATAC.3' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.3
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/SHARE-seq/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'BloodDev'){
    if('SeuratObj_ATAC.4' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.4
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/BloodDev/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'Kidney'){
    if('SeuratObj_ATAC.5' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.5
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Kidney/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'Brain'){
    if('SeuratObj_ATAC.6' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.6
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/Brain/SeuratObj/SeuratObj_ATAC.RDS')
    }
  } else if(dataset == 'ISSAAC-seq'){
    if('SeuratObj_ATAC.12' %in% ls){
      SeuratObj_ATAC <- SeuratObj_ATAC.12
    } else{
      SeuratObj_ATAC <- readRDS('~/txmdata/scATAC_seq/AtacAnnoR/data/ISSAAC-seq/SeuratObj/SeuratObj_ATAC.RDS')
    }
  }
  return(SeuratObj_ATAC)
}



