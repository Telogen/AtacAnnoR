################### GetMergedLabels ##################
# get_merged_labels_Kidney
# get_merged_labels_Zhu
# get_merged_labels_Wilk
# get_merged_labels_Stephenson
# get_merged_labels_Hao
# get_merged_labels_Monaco


get_merged_labels_dataset_ISSAACseq <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels ==      "R0 Ex-L2/3 IT")] = "R0"
    labels[which(labels ==  "R1 Ex-L2/3 IT Act")] = "R1"
    labels[which(labels ==         "R10 Ex-L6b")] = "R10"
    labels[which(labels ==   "R11 Ex-PIR Ndst4")] = "R11"
    labels[which(labels ==        "R13 In-Drd2")] = "R13"
    labels[which(labels ==        "R14 In-Hap1")] = "R14"
    labels[which(labels ==       "R15 In-Pvalb")] = "R15"
    labels[which(labels ==         "R16 In-Sst")] = "R16"
    labels[which(labels ==        "R17 In-Tac1")] = "R17"
    labels[which(labels ==   "R18 In-Vip/Lamp5")] = "R18"
    labels[which(labels ==          "R19 Astro")] = "R19"
    labels[which(labels ==        "R2 Ex-L4 IT")] = "R2"
    labels[which(labels ==            "R20 OPC")] = "R20"
    labels[which(labels ==          "R21 Oligo")] = "R21"
    labels[which(labels ==           "R22 VLMC")] = "R22"
    labels[which(labels ==        "R3 Ex-L5 IT")] = "R3"
    labels[which(labels ==        "R4 Ex-L5 NP")] = "R4"
    labels[which(labels == "R5 Ex-L5 NP Cxcl14")] = "R5"
    labels[which(labels ==        "R6 Ex-L5-PT")] = "R6"
    labels[which(labels ==        "R7 Ex-L6 CT")] = "R7"
    labels[which(labels ==   "R8 Ex-L6 IT Bmp3")] = "R8"
    labels[which(labels ==  "R9 Ex-L6 IT Oprk1")] = "R9"
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}



get_merged_labels_dataset_SHAREseq <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels == "ahighCD34+ bulge"         )] = 'HCD34B'
    labels[which(labels == "alowCD34+ bulge"          )] = 'LCD34B'
    labels[which(labels == "Basal"                    )] = 'Bas'
    labels[which(labels == "Dermal Fibroblast"        )] = 'DF'
    labels[which(labels == "Dermal Papilla"           )] = 'DP'
    labels[which(labels == "Dermal Sheath"            )] = 'DS'
    labels[which(labels == "Endothelial"              )] = 'Endo'
    labels[which(labels == "Granular"                 )] = 'Gran'
    labels[which(labels == "Hair Shaft-cuticle.cortex")] = 'HSCC'
    labels[which(labels == "Infundibulum"             )] = 'Infu'
    labels[which(labels == "IRS"                      )] = 'IRS'
    labels[which(labels == "Isthmus"                  )] = 'Isth'
    labels[which(labels == "K6+ Bulge Companion Layer")] = 'KBCL'
    labels[which(labels == "Macrophage DC"            )] = 'MDC'
    labels[which(labels == "Medulla"                  )] = 'Medu'
    labels[which(labels == "Melanocyte"               )] = 'Mela'
    labels[which(labels == "ORS"                      )] = 'ORS'
    labels[which(labels == "Schwann Cell"             )] = 'SC'
    labels[which(labels == "Sebaceous Gland"          )] = 'SG'
    labels[which(labels == "Spinous"                  )] = 'Spin'
    labels[which(labels == "TAC-1"                    )] = 'TAC1'
    labels[which(labels == "TAC-2"                    )] = 'TAC2'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}



get_merged_labels_Kidney <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels == 'PCT')] ='PT'
    labels[which(labels == 'PST')] ='PT'
    labels[which(labels == 'DCT1')]='DCT'
    labels[which(labels == 'DCT2')]='DCT'
    labels[which(labels == 'MES')] ='MES_FIB'
    labels[which(labels == 'FIB')] ='MES_FIB'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}


get_merged_labels_Zhu <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD16 Mono","CD14 Mono","Monocytes"))] = 'Mono'
    labels[which(labels %in% c("cDC", "pDC","DCs"))] = 'DC'
    labels[which(labels %in% c("CD4 Naive", "CD8 Naive", "Naive T cells"))] = 'NaiveT'
    labels[which(labels %in% c("CD4 TCM", "Treg", "CD4 TEM", "Activated CD4 T cells"))] = 'CD4T'
    labels[which(labels %in% c("CD8 TEM_2", "CD8 TEM_1", "Cytotoxic CD8 T cells"))] = 'CD8T'
    labels[which(labels %in% c("NK", "NKs", "XCL+ NKs"))] = 'ILC'
    labels[which(labels %in% c("Memory B", "Memory B cells","Intermediate B"))] = 'MemB'
    labels[which(labels %in% c("Naive B", "Naive B cells"))] = 'NaiveB'
    labels[which(labels %in% c("Plasma", "Cycling Plasma"))] = 'Plasma'
    labels[which(labels %in% c("Cycling T cells"))] = 'CycT'
    labels[which(labels %in% c("Megakaryocytes"))] = 'Mega'
    labels[which(labels %in% c("Stem cells"))] = 'HSPC'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}



get_merged_labels_Wilk <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD16 Mono", "CD16 Monocyte"))] = 'CD16Mono'
    labels[which(labels %in% c("CD14 Mono", "CD14 Monocyte"))] = 'CD14Mono'
    labels[which(labels %in% c("cDC", "DC"))] = 'cDC'
    labels[which(labels %in% c("NK"))] = 'ILC'
    labels[which(labels %in% c("CD4 Naive", "CD8 Naive", "CD4n T"))] = 'NaiveT'
    labels[which(labels %in% c("CD4 TCM", "Treg", "CD4 TEM", "CD4m T","CD4 T"))] = 'CD4T'
    labels[which(labels %in% c("CD8m T", "CD8 TEM_2", "CD8 TEM_1","MAIT","CD8eff T"))] = 'CD8T'
    labels[which(labels %in% c("gdT", "gd T"))] = 'gdT'
    labels[which(labels %in% c("Intermediate B", "Memory B", "Naive B","B"))] = 'B'
    labels[which(labels %in% c("Plasmablast", "Plasma"))] = 'Plasma'
    labels[which(labels %in% c("SC & Eosinophil", "HSPC"))] = 'HSPC'
    labels[which(labels %in% c("Granulocyte"))] = 'Granu'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}



get_merged_labels_Stephenson <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD16 Mono", "CD16.mono"))] = 'CD16Mono'
    labels[which(labels %in% c("CD14 Mono", "CD14.mono"))] = 'CD14Mono'
    labels[which(labels %in% c("cDC", "DC"))] = 'cDC'
    labels[which(labels %in% c("CD4 Naive", "CD4.Naive"))] = 'CD4NT'
    labels[which(labels %in% c("CD8 Naive", "CD8.Naive"))] = 'CD8NT'
    labels[which(labels %in% c("CD4 TCM", "CD4 TEM","CD4.CM", "CD4.IL22", "CD4.Th",
                               "CD4.EM", "CD4.Tfh"))] = 'CD4T'
    labels[which(labels %in% c("CD8.TE", "CD8.EM", "CD8 TEM_2", "CD8 TEM_1"))] = 'CD8T'
    labels[which(labels %in% c("Intermediate B", "Memory B","B_non-switched_memory",
                               "B_switched_memory","B_exhausted"))] = 'MemB'
    labels[which(labels %in% c("B_naive", "Naive B","B_immature"))] = 'NaiveB'
    labels[which(labels %in% c("HSC", "HSPC"))] = 'HSPC'
    labels[which(labels %in% c("NK", "ILC"))] = 'ILC'
    labels[which(labels %in% c("Lymph.prolif"))] = 'Prolif'
    labels[which(labels %in% c("Platelets"))] = 'Platelet'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}





get_merged_labels_Hao <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD16 Mono"))] = 'CD16Mono'
    labels[which(labels %in% c("CD14 Mono"))] = 'CD14Mono'
    labels[which(labels %in% c("CD4 Naive"))] = 'CD4NT'
    labels[which(labels %in% c("CD8 Naive"))] = 'CD8NT'
    labels[which(labels %in% c("CD4 CTL"))] = 'CD4CTL'
    labels[which(labels %in% c("CD8 TEM","CD8 TCM", "CD8 TEM_2", "CD8 TEM_1"))] = 'CD8T'
    labels[which(labels %in% c("CD4 TEM","CD4 TCM","CD4 CTL"))] = 'CD4T'
    labels[which(labels %in% c("Intermediate B", "B intermediate"))] = 'InterB'
    labels[which(labels %in% c("Memory B", "B memory"))] = 'MemB'
    labels[which(labels %in% c("Naive B", "B naive"))] = 'NaiveB'
    labels[which(labels %in% c("Plasmablast", "Plasma"))] = 'Plasma'
    labels[which(labels %in% c("NK", "ILC"))] = 'ILC'
    labels[which(labels %in% c("pDC","ASDC"))] = 'pDC'
    labels[which(labels %in% c("Proliferating"))] = 'Prolif'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}




get_merged_labels_Monaco <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD16 Mono", "NC_mono"))] = 'CD16Mono'
    labels[which(labels %in% c("CD14 Mono", "C_mono"))] = 'CD14Mono'
    labels[which(labels %in% c("I_mono"))] = 'InterMono'
    labels[which(labels %in% c("NK"))] = 'ILC'
    labels[which(labels %in% c("cDC", "mDC"))] = 'cDC'
    labels[which(labels %in% c("CD4 Naive", "CD4_naive"))] = 'CD4NT'
    labels[which(labels %in% c("CD8 Naive", "CD8_naive"))] = 'CD8NT'
    labels[which(labels %in% c("CD4 TCM", "CD4 TEM", "CD4_TE","TFH", "Th1",
                               "Th1.Th17", "Th17", "Th2","Th1/Th17"))] = 'CD4T'
    labels[which(labels %in% c( "CD8 TEM_2", "CD8 TEM_1","CD8_CM", "CD8_EM", "CD8_TE"))] = 'CD8T'
    labels[which(labels %in% c("gdT", "VD2-","VD2+","VD2_gdT","nVD2_gdT"))] = 'gdT'
    labels[which(labels %in% c("Intermediate B", "Memory B","B_NSM","B_Ex", "B_SM"))] = 'MemB'
    labels[which(labels %in% c("Naive B", "B_naive"))] = 'NaiveB'
    labels[which(labels %in% c("Plasmablasts", "Plasma"))] = 'Plasma'
    labels[which(labels %in% c("Progenitor", "HSPC"))] = 'HSPC'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}




get_Zhu_score <- function(celltype){
  score = c(1,1,2,2,1,2,2,2,1,1,1,1,1,1,1)
  names(score) <- c('Naive T cells','Cytotoxic CD8 T cells','MAIT','NKs','Activated CD4 T cells',
                    'Naive B cells','Memory B cells','XCL+ NKs','Cycling T cells','Monocytes',
                    'DCs','Megakaryocytes','Plasma','Cycling Plasma','Stem cells')
  return(score[celltype])
}

get_Wilk_score <- function(celltype){
  score = c(1,2,2,1,1,1,2,1,2,1,2,2,1,1,1,1,1)
  names(score) <- c('NK','CD8m T','CD14 Monocyte','B','CD4n T','CD4m T','CD16 Monocyte',
                    'gd T','DC','RBC','pDC','CD8eff T','Plasmablast','SC & Eosinophil',
                    'Platelet','Granulocyte','CD4 T')
  return(score[celltype])
}


get_Stephenson_score <- function(celltype){
  score = c(2,1,2,2,1,2,2,2,2,1,2,2,2,1,3,2,1,2,3,3,1,1,1,2,3,2,1,2)
  names(score) <- c('NK','CD4.Naive','CD14.mono','CD4.CM','CD8.Naive','CD4.IL22','CD8.TE',
                    'CD8.EM','B_naive','gdT','MAIT','CD16.mono','DC','Platelets',
                    'B_switched_memory','CD4.Th','NKT','pDC','B_immature','B_non-switched_memory',
                    'Plasma','Lymph.prolif','RBC','ILC','B_exhausted','CD4.EM','HSC','Treg')
  return(score[celltype])
}

get_Hao_score <- function(celltype){
  score = c(2,2,2,2,2,2,2,2,2,1,3,2,2,2,2,3,1,2,2,1,1,1,1,2,1,3)
  names(score) <- c('CD14 Mono','NK','CD4 Naive','CD4 TCM','CD8 TEM','CD8 Naive','B naive',
                    'CD16 Mono','CD4 TEM','gdT','B memory','CD8 TCM','MAIT','cDC',
                    'Treg','B intermediate','Platelet','CD4 CTL','pDC','Proliferating',
                    'Plasmablast','dnT','HSPC','ILC','Eryth','ASDC')
  return(score[celltype])
}

get_Monaco_score <- function(celltype){
  score = c(3,2,3,3,1,2,2,2,2,2,2,2,2,2,2,1,1,2,2,1,1,2,2,2,2,2,2,2,2)
  names(score) <- c('B_Ex','B_naive','B_NSM','B_SM','Basophils','C_mono','CD4_naive',
                    'CD8_CM','CD8_EM','CD8_naive','CD8_TE','I_mono','MAIT','mDC',
                    'NC_mono','Neutrophils','NK','nVD2_gdT','pDC','Plasmablasts',
                    'Progenitor','TFH','Th1','Th1.Th17','Th17','Th2','Treg',
                    'VD2_gdT','CD4_TE')
  return(score[celltype])
}



get_merged_labels_dataset0 <- function(df){
  if(class(df) != 'data.frame'){
    df_new <- data.frame(df = df)
  } else{
    df_new <- df
  }
  df_new <- apply(df_new,2,function(labels){
    labels <- as.character(labels)
    labels[which(labels %in% c("CD14 Mono"))] = 'CD14Mono'
    labels[which(labels %in% c("CD16 Mono"))] = 'CD16Mono'
    labels[which(labels %in% c("CD4 Naive"))] = 'CD4NT'
    labels[which(labels %in% c("CD8 Naive"))] = 'CD8NT'
    labels[which(labels %in% c("CD4 TCM"))] = 'CD4TCM'
    labels[which(labels %in% c("CD4 TEM"))] = 'CD4TEM'
    labels[which(labels %in% c("NK"))] = 'ILC'
    labels[which(labels %in% c("Naive B"))] = 'NaiveB'
    labels[which(labels %in% c("CD8 TEM_2", "CD8 TEM_1"))] = 'CD8T'
    labels[which(labels %in% c("Memory B"))] = 'MemB'
    labels[which(labels %in% c("Intermediate B"))] = 'InterB'
    return(labels)
  }) %>% as.data.frame()
  rownames(df_new) <- rownames(df)
  if(class(df) != 'data.frame'){
    return(df_new[,1])
  } else{
    return(df_new)
  }
}

get_dataset0_score <- function(celltype){
  score = c(2,2,2,2,2,1,2,2,3,2,2,2,2,1,2,2,2,1,1)
  names(score) <- c('CD14 Mono','CD4 Naive','CD8 Naive','CD4 TCM','CD16 Mono','NK',
                    'Naive B','CD8 TEM_2','Intermediate B','CD8 TEM_1','CD4 TEM',
                    'cDC','Treg','gdT','Memory B','MAIT','pDC','HSPC','Plasma')
  return(score[celltype])
}

