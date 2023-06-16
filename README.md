# AtacAnnoR: A Reference-Based Annotation Tool for Single Cell ATAC-seq Data

## Overview

AtacAnnoR is a novel scATAC-seq cell type annotation method using scRNA-seq data as the reference. 

<img src="https://github.com/Telogen/AtacAnnoR/blob/main/figures/Fig1.png" width="800">

- AtacAnnoR performs **two rounds** of annotation, which annotate scATAC-seq cells at gene-level and genome-wide-level, respectively.

- In the first round **gene-level annotation**, AtacAnnoR compares gene activity profiles derived from scATAC-seq data with reference gene expression profiles, assigning a subset of query cells with reference cell type labels. 

- In the second round **genome-wide-level annotation**, AtacAnnoR predicts the labels of the remaining query cells using a meta-program matrix derived from genome-wide ATAC-peaks. 



## Installation

```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("zdebruine/RcppML")        # v0.5.4 or higher
devtools::install_github("TianLab-Bioinfo/AtacAnnoR",force = T)
```


## Usage


### Run AtacAnnoR in a line with default parameters

```
pred <- RunAtacAnnoR(ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                     ref_celltype = SeuratObj_RNA$true, 
                     query_gene_activity = SeuratObj_ATAC[['ACTIVITY']]@counts, 
                     query_peak_counts = SeuratObj_ATAC[['ATAC']]@counts) 
```

### Run AtacAnnoR step by step

- You can run AtacAnnoR step by step to see how AtacAnnoR works and modify parameters to annotate scATAC-seq cells better. See [this tutorial](https://telogen.github.io/AtacAnnoR/Run_AtacAnnoR_step_by_step_v2.html).


### AtacAnnoR “Combine and Discard” strategy

- When multiple references' annotation results are available, AtacAnnoR applies a “Combine and Discard” strategy to discard some cells and further improve the annotation accuracy. See [this tutorial](https://telogen.github.io/AtacAnnoR/Combine_and_Discard.html).



### Beta function: AtacAnnoR label corrector

- AtacAnnoR also provides a more aggressive strategy to use the multiple references' annotation results to correct the original labels and improve accuracy. See [this tutorial](https://telogen.github.io/AtacAnnoR/AtacAnnoR_label_corrector.html).




### Run AtacAnnoR in other scATAC-seq analysis pipelines


- Run AtacAnnoR in [Signac](https://stuartlab.org/signac)

```
SeuratObj_ATAC <- RunAtacAnnoR_Signac(query_SeuratObj = SeuratObj_ATAC,
                                      query_ga_assay = 'ACTIVITY',
                                      query_peak_assay = 'ATAC',
                                      ref_SeuratObj = SeuratObj_RNA,
                                      ref_assay = 'RNA',
                                      ref_ident = 'true')
```



- Run AtacAnnoR in [ArchR](https://www.archrproject.com/bookdown)

```
query_ArchRproj <- RunAtacAnnoR_ArchR(query_ArchRproj = query_ArchRproj,
                                      ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                                      ref_celltype = SeuratObj_RNA$true)
```


- Run AtacAnnoR in [SnapATAC](https://github.com/r3fang/SnapATAC)


```
query_snapObj <- RunAtacAnnoR_SnapATAC(query_snapObj = query_snapObj,
                                       ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                                       ref_celltype = SeuratObj_RNA$true)
```


- Run AtacAnnoR in [Cicero](https://cole-trapnell-lab.github.io/cicero-release/)

```
query_cds <- RunAtacAnnoR_Cicero(query_cds = query_cds,
                                 ref_mtx = SeuratObj_RNA[['RNA']]@counts, 
                                 ref_celltype = SeuratObj_RNA$true, 
                                 query_gene_activity = query_gene_activity)
```



## Contact

Lejin Tian: ljtian20@fudan.edu.cn


