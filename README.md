# <Integrative multi-omics and machine learning approaches uncover a novel metabolic-related signature associated with cancer-associated fibroblasts in gastric cancer development>

> **Paper title: Integrative multi-omics and machine learning approaches uncover a novel metabolic-related signature associated with cancer-associated fibroblasts in gastric cancer development** 



## Table of Contents

- [About](#about)  
- [Repository Structure](#repository-structure)  
- [Prerequisites](#prerequisites)   
- [Step by step](#step-by-step)



## About

This pipeline performs:

1. **Single-cell RNA-seq**  
   - Identification of CAF clusters with high expression of metabolic signatures.

2. **Bulk analysis**  
   - Integration of multiple GEO series and TCGA-STAD data.

3. **Machine Learning**  
   - Leave-one-out cross-validation of Cox-based and survival models to establish MRS.

4. **Evaluation of MRS**  
   - Generation of forest plots and nomograms for clinical variables and risk scores.

5. **Somatic mutation**  
   - Extraction of MAF files to stratify patients by MRS.


High-quality, reproducible analyses with minimal manual downloads.



--- 

##  Repository Structure

```text
├── scripts/
│   ├── TCGA-STAD-Lisapaper-revise.R
│   ├── ML_input_preparation.R
│   ├── Machine_gastric_loocv.R
│   ├── gse183904_gastric_SINGLECELL.R
│   └── forestplot_nomogram_clinical.R
|   └── MAF_file_gastric.R
└── README.md
```
---
##  Prerequisites
R ≥ 4.1.0

R packages (tested versions in parentheses):

- `circlize` (0.4.16)  
- `clusterProfiler` (4.14.3)  
- `dbplyr` (2.5.0)  
- `DESeq2` (1.46.0)  
- `dplyr` (1.1.4)  
- `fgsea` (1.32.0)  
- `GenomeInfoDb` (1.42.0)  
- `ggplot2` (3.5.1)  
- `GSVA` (2.0.1)  
- `igraph` (2.1.1)  
- `IOBR` (0.99.0)  
- `maftools` (2.22.0)  
- `magick` (2.8.5)  
- `Matrix` (1.7-1)  
- `patchwork` (1.3.0)  
- `pheatmap` (1.0.12)  
- `plyr` (1.8.9)  
- `R.utils` (2.12.3)  
- `RColorBrewer` (1.1-3)  
- `sessioninfo` (1.2.2)  
- `Seurat` (5.1.0)  
- `SummarizedExperiment` (1.36.0)  
- `survival` (3.7-0)  
- `survminer` (0.5.0)  
- `TCGAbiolinks` (2.34.0)  
- `tidyHeatmap` (1.8.1)  
- `tidyverse` (2.0.0)  
- `timeROC` (0.4)  

---
##  Step by step

### Single-cell analysis
- **Script:** `gse183904_gastric_SINGLECELL.R`  
  _Generates Figures 1 and 2_

### Bulk & DE analysis
- **Script:** `ML_input_preparation.R`  
  _Prepares inputs for bulk RNA-seq and microarray analysis_

### Machine learning & MRS construction
- **Script:** `TCGA-STAD-Lisapaper-revise.R`  
  _Generates Figure 3A_
- **Script:** `Machine_gastric_loocv.R`  
  _Generates Figures 3B, 3C, 3D_

### Evaluation of MRS performance
- **Script:** `Machine_gastric_loocv.R`  
  _Generates Figures 4A–4F_
- **Script:** `forestplot_nomogram_clinical.R`  
  _Generates Figures 4G–4J and 5_
- **Script:** `MAF_file_gastric.R`  
  _Generates Figure 6_
---
