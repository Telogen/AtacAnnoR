# AtacAnnoR

## Overview

AtacAnnoR is a method to annotate scATAC-seq cell type from well-annotated scRNA-seq data.

<img src="https://github.com/Telogen/AtacAnnoR/blob/main/figures/pipeline.png" width="600">

## Installation

- Required packages:
```R
install.packages('Matrix')
install.packages('Seuart')
install.packages('pcaPP')
install.packages('pbmcapply')
devtools::install_github("zdebruine/RcppML")        # v0.5.4 or higher
devtools::install_github("Telogen/AtacAnnoR")
```

- Recommended packages: (will not affect the use of basic functions)

```R
install.packages('ComplexHeatmap')
install.packages('circlize')
install.packages('grid')
install.packages('ggplot2')
install.packages('cowplot')
install.packages('aplot')
```

## Tutorial

### Get the input for AtacAnnoR

- We recommend to use [Seurat](https://satijalab.org/seurat/index.html) to process scRNA-seq data and [Signac](https://satijalab.org/signac/) to process scATAC-seq data. You can get the input for AtacAnnoR following [this tutorial](https://github.com/Telogen/AtacAnnoR/blob/main/tutorial/get_scMAGICatac_input.ipynb).

### A quick start to AtacAnnoR

- Run AtacAnnoR in a line with default parameters, see [this tutorial](https://github.com/Telogen/AtacAnnoR/blob/main/tutorial/quick_start_to_scMAGICatac.ipynb).


### Run AtacAnnoR step by step

- You can also run AtacAnnoR step by step to see how AtacAnnoR works and modify parameters to annotate your scATAC-seq cells better, see [this tutorial](https://github.com/Telogen/AtacAnnoR/blob/main/tutorial/run_scMAGICatac_step_by_step.ipynb).


### Run AtacAnnoR in Signac, SnapATAC, ArchR and Cicero pipeline

- See [this tutorial](https://github.com/Telogen/AtacAnnoR/blob/main/tutorial/run_scMAGICatac_in_other_pipelines.ipynb).

## Citation

- coming soon

## Contact

- Lejin Tian: ljtian20@fudan.edu.cn


