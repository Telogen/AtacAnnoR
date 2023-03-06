# AtacAnnoR

## Overview

AtacAnnoR: A Reference-Based Annotation Tool for Single Cell ATAC-seq Data

<img src="https://github.com/Telogen/AtacAnnoR/blob/main/figures/fig1A.png" width="800">

## Installation

```R
devtools::install_github("zdebruine/RcppML")        # v0.5.4 or higher
devtools::install_github("TianLab-Bioinfo/AtacAnnoR")
```


## Usage

```
pred <- RunAtacAnnoR(ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                     ref_celltype = SeuratObj_RNA$true, 
                     ref_type = "sc",
                     query_gene_activity = SeuratObj_ATAC[['ACTIVITY']]@counts, 
                     query_peak_counts = SeuratObj_ATAC[['ATAC']]@counts, 
                     query_nmf_embedding = NULL,
                     threads = 10, 
                     verbose = TRUE, 
                     simple_output = TRUE) 

get_benchmark(SeuratObj_ATAC$true,pred)
#          accuracy    average_recall average_precision          macro_f1 
#         0.9132731         0.8897061         0.8968420         0.8904343 
```

- Run AtacAnnoR in Signac pipeline

```
SeuratObj_ATAC <- RunAtacAnnoR_Signac(ref_SeuratObj = SeuratObj_RNA,
                                      ref_assay = 'RNA',
                                      ref_ident = 'true',
                                      ref_type = "sc",
                                      query_SeuratObj = SeuratObj_ATAC,
                                      query_ga_assay = 'ACTIVITY',
                                      query_peak_assay = 'ATAC',
                                      threads = 10, 
                                      verbose = TRUE)
```

- Run AtacAnnoR in SnapATAC pipeline

```
query_snapObj <- RunAtacAnnoR_SnapATAC(ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                                       ref_celltype = SeuratObj_RNA$true, 
                                       ref_type = "sc",
                                       query_snapObj = query_snapObj,
                                       threads = 10, 
                                       verbose = TRUE)
```


- Run AtacAnnoR in ArchR pipeline

```
query_ArchRproj <- RunAtacAnnoR_ArchR(ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                                      ref_celltype = SeuratObj_RNA$true, 
                                      ref_type = "sc",
                                      query_ArchRproj = query_ArchRproj,
                                      threads = 10, 
                                      verbose = TRUE)
```



## Contact

- Lejin Tian: ljtian20@fudan.edu.cn


