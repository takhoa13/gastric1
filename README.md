# <Integrative multi-omics and machine learning approaches uncover a novel metabolic-related signature associated with cancer-associated fibroblasts in gastric cancer development>

> **One-sentence summary:** <Short description of what this project does>

---

## Table of Contents

- [About](#about)  
- [Repository Structure](#repository-structure)  
- [Prerequisites](#prerequisites)   
- [Step by step](#scripts)  
- [Contributing](#contributing)  
- [Contact](#contact)  

---

## About

This pipeline performs:

1. **Bulk RNA-seq (TCGA-STAD)**  
   - Differential expression analysis & Cox regression.
2. **Microarray (GEO)**  
   - Integration of multiple GEO series, feature selection.
3. **Machine Learning**  
   - Leave-one-out cross-validation of Cox-based and tree-based survival models.
4. **Single-cell RNA-seq**  
   - QC, integration (Seurat + Harmony), clustering.
5. **Visualization**  
   - Forest plots and nomograms for clinical variables and risk scores.

High-quality, reproducible analyses with minimal manual downloads.


---

## 🗂️ Repository Structure

```text
├── scripts/
│   ├── TCGA-STAD-Lisapaper-revise.R
│   ├── ML_input_preparation.R
│   ├── Machine_gastric_loocv.R
│   ├── gse183904_gastric_SINGLECELL.R
│   └── forestplot_nomogram_clinical.R
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
**Single-cell analysis**
`gse183904_gastric_SINGLECELL.R` refers to Figure 1 and Figure 2
**Bulk & DE analysis**
`ML_input_preparation.R` refers to the preparation steps for bulk analysis
**Machine learning and establishing MRS**
`TCGA-STAD-Lisapaper-revise.R` refers to Figure 3A
`Machine_gastric_loocv.R` refers to Figure 3B, 3C, 3D
**Evaluation of MRS performance**
`Machine_gastric_loocv.R` refers to Figure 4A-F
`forestplot_nomogram_clinical.R` refers to Figure 4G-J and Figure 5
`MAF_file_gastric.R` refers to Figure 6
