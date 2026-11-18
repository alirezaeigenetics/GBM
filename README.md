# Glioblastoma (GBM) Transcriptomic Meta-Analysis & Machine Learning Pipeline

This repository contains the complete, reproducible R and Python computational pipeline for discovering a universally detected, non-invasive 4-microRNA diagnostic candidate panel (`hsa-miR-137-3p`, `hsa-miR-769-3p`, `hsa-miR-485-3p`, and `hsa-miR-765`) for Glioblastoma (GBM).

To overcome single-cohort biases, the pipeline integrates six independent microarray datasets from NCBI GEO and one RNA-sequencing cohort from TCGA. It automates data acquisition, quality control, cross-platform normalization, differential expression analysis, Fisher's meta-analysis, survival modeling (Exo-GBM Prognostic Index, EGPI), topologically-weighted machine learning, and Explainable AI (SHAP).

---

## Repository Structure

```text
GBM/
├── Main.Rproj                    # Main RStudio Project file (OPEN THIS FIRST)
├── README.md                     # Project documentation
├── LICENSE                       # Repository license
├── Merge.py                      # Python script to merge all *_ML_Ready.csv files
├── ML.py                         # Python script for Topologically-Weighted ML, Ablation Study & SHAP
├── GSE61710/                     # GEO microarray dataset & analysis script
├── GSE65626/                     # GEO microarray dataset & analysis script
├── GSE90603/                     # GEO microarray dataset & analysis script
├── GSE91014/                     # GEO microarray dataset & analysis script
├── GSE93850/                     # GEO microarray dataset & analysis script
├── GSE103228/                    # GEO microarray dataset & analysis script
├── TCGA_RNAseq/                  # TCGA-GBM RNA-seq, DESeq2 & Survival (EGPI) script
└── Meta Analysis/                # Cross-cohort Venn, UpSet & Fisher's Meta-Analysis script
```

---

## Prerequisites

* **R (v4.5.1 or higher)** & **RStudio**: Required R and Bioconductor packages (`BiocManager`, `limma`, `DESeq2`, `oligo`, `TCGAbiolinks`, `GEOquery`, `survival`, `survminer`, `ggplot2`, `ggrepel`, `ggVennDiagram`, `UpSetR`, etc.) are automatically checked and installed at the top of each R script.
* **Python (v3.10 or higher)**: Install the required libraries before running the Python scripts:
  ```bash
  pip install pandas numpy scikit-learn shap matplotlib scipy
  ```

---

## Step-by-Step Execution Workflow

### Step 1: Individual Dataset Processing (R)
> **Important:** To ensure all relative file paths resolve correctly, you **must** double-click and open `Main.Rproj` in RStudio first. Do not run the scripts standalone without loading the R project.

1. Open `Main.Rproj` in RStudio.
2. Open and run the R script inside each of the seven dataset folders (`GSE61710/`, `GSE65626/`, `GSE90603/`, `GSE91014/`, `GSE93850/`, `GSE103228/`, and `TCGA_RNAseq/`).
3. Each script will automatically:
   * Download raw data and clinical metadata from GEO or TCGA-GDC.
   * Perform platform-specific normalization (RMA, quantile normalization, or `rlog`) and differential expression analysis (`limma` or `DESeq2`).
   * Save high-resolution vector visualizations (`.svg`), including Boxplot, PCA, UMAP, Mean-Variance, Mean-Difference, Quantile-Quantile plot, Volcano plot, Adjusted p-value histogram, Venn diagram, Heatmap, Density plot, and Scatter plot.
   * Export a Z-score normalized `*_ML_Ready.csv` file inside each dataset folder.
   * *(TCGA only)* Run Kaplan-Meier survival analysis with MaxStat optimal cutpoints and compute the multivariate **Exo-GBM Prognostic Index (EGPI)**.

### Step 2: Cross-Cohort Meta-Analysis (R)
1. While still inside the `Main.Rproj` RStudio session, open and run the R script inside the `Meta Analysis/` folder.
2. This script reads the differential expression results from all seven dataset folders, standardizes miRNA nomenclature, and generates:
   * Overlap visualizations (`MetaAnalysis_All_VennDiagram` and `MetaAnalysis_All_UpSetPlot` in `.png` and `.svg`).
   * A list of universally detected miRNAs (`Conserved_All_miRNAs.csv`).
   * Fisher's Combined Probability Test results (`MetaAnalysis_Fisher_Results.csv`) and a labeled meta-volcano plot (`MetaAnalysis_Custom_Volcano_Labeled`).

### Step 3: Merge Datasets for Machine Learning (Python)
Run `Merge.py` from the root directory:
```bash
python Merge.py
```
This script recursively scans all dataset subfolders for `*_ML_Ready.csv` files, performs an inner join to retain only universally detected (conserved) miRNAs across all seven cohorts, shuffles the rows, and outputs `MASTER_GBM_DATASET.csv` in the root folder.

### Step 4: Machine Learning, Ablation Study & Explainable AI (Python)
Run `ML.py` from the root directory:
```bash
python ML.py
```
This script loads `MASTER_GBM_DATASET.csv` and executes:
* **10-Iteration Balanced Random Undersampling** with **5-Fold Stratified Cross-Validation**.
* **Ablation Study**: Compares 5 feature selection configurations (`RF_Alone`, `LASSO_Alone`, `SVM_Alone`, `Consensus_Unweighted`, and `Consensus_Weighted` using network degree centrality) and prints AUC, Accuracy, Sensitivity, and Specificity to the console.
* **Explainable AI (SHAP)**: Saves `SHAP_Summary_Plot.png` and individual `SHAP_Dependence_*.png` plots for the 4-miRNA panel.
* **ROC Curve**: Saves `High_Def_ROC_Curve.png` (600 DPI).

---

## Notes on Network Access
Automated downloads from NCBI GEO and the NCI GDC portal (`GEOquery` and `TCGAbiolinks`) require a stable internet connection. If you encounter connection timeouts or regional network restrictions, you can route R through a local proxy or manually download the raw archives into the corresponding dataset folder before running the scripts.

## License
This project is licensed under the terms specified in the `LICENSE` file.
