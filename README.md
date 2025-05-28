# <Integrative multi-omics and machine learning approaches uncover a novel metabolic-related signature associated with cancer-associated fibroblasts in gastric cancer development>

> **One-sentence summary:** <Short description of what this project does>

---

## Table of Contents

- [About](#about)  
- [Repository Structure](#repository-structure)  
- [Prerequisites](#prerequisites)   
- [Scripts](#scripts)  
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
