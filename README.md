# Glioblastoma (GBM) Transcriptomic Meta-Analysis & Machine Learning Pipeline

This repository contains reproducible R scripts for the automated download, quality control, normalization, and differential expression analysis of multiple Glioblastoma (GBM) microarray and RNA-seq datasets from NCBI GEO and TCGA. 

The pipeline processes raw transcriptomic data into normalized expression matrices and identifies significant differentially expressed genes (DEGs). It also outputs standardized, Z-score normalized datasets formatted for downstream Machine Learning (ML) integration.

## Repository Structure

```text
GBM/
├── Main.Rproj                     # Main RStudio Project file (Open this first)
├── README.md                     # Project documentation
├── Combine_ML_Datasets.R         # Master script to merge all ML_ready.csv files (Example)
├── GSE61710/                     # Dataset-specific directory
│   ├── GSE61710.R                # Primary analysis script for GSE61710
│   └── (Outputs generate here)
└── GSE90603/                     # Dataset-specific directory
    ├── GSE90603.R                # Primary analysis script for GSE90603
    └── (Outputs generate here)
