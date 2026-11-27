# Glioblastoma (GBM) Transcriptomic Meta-Analysis & Machine Learning Pipeline

This repository contains the complete, reproducible R and Python computational pipeline for discovering and validating a universally detected, tissue-derived 4-microRNA diagnostic candidate panel (`hsa-miR-137-3p`, `hsa-miR-769-3p`, `hsa-miR-485-3p`, and `hsa-miR-765`) with extracellular vesicle (EV) translational potential for Glioblastoma (GBM).

To overcome single-cohort biases, the pipeline integrates six independent microarray datasets from NCBI GEO and one RNA-sequencing cohort from TCGA. It automates data acquisition, quality control, cross-platform normalization, differential expression analysis, Fisher's meta-analysis, survival modeling (GBM miRNA Prognostic Index, GMPI) with clinical comparator benchmarking, topologically-weighted machine learning with paired statistical significance testing, 7-fold Leave-One-Dataset-Out (LODO) cross-validation, and Explainable AI (SHAP).

---

## Repository Structure

```text
GBM/
├── Main.Rproj                     # Main RStudio Project file (OPEN THIS FIRST)
├── README.md                      # Project documentation
├── LICENSE                        # Repository license
├── ML_pipeline.py                 # Unified Python pipeline (Merge, CV, Ablation, LODO & SHAP)
├── Ablation_Study_Results.csv     # Auto-generated 10x5-fold CV ablation table (95% CIs & p-values)
├── LODO_Validation_Results.csv    # Auto-generated 7-fold & pooled LODO validation summary
├── ML_Outputs/                    # Auto-generated folder for high-res (.png & .svg) figures
├── GSE61710/                      # GEO microarray dataset & analysis script (peritumoral excluded for ML)
├── GSE65626/                      # GEO microarray dataset & analysis script
├── GSE90603/                      # GEO microarray dataset & analysis script
├── GSE91014/                      # GEO microarray dataset & analysis script
├── GSE93850/                      # GEO microarray dataset & analysis script
├── GSE103228/                     # GEO microarray dataset & analysis script
├── TCGA_RNAseq/                   # TCGA-GBM RNA-seq, DESeq2, Survival (GMPI) & Clinical Comparator script
└── Meta Analysis/                 # Cross-cohort Venn, UpSet & Fisher's Meta-Analysis script
```

---

## Prerequisites

- **R (v4.5.1 or higher)** & **RStudio**: Required R and Bioconductor packages (`BiocManager`, `limma`, `DESeq2`, `oligo`, `TCGAbiolinks`, `GEOquery`, `survival`, `survminer`, `timeROC`, `dcurves`, `ggplot2`, `ggrepel`, `ggVennDiagram`, `UpSetR`, etc.) are automatically checked and installed at the top of each R script.
- **Python (v3.10 or higher)**: Install the required libraries before running the Python pipeline:
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
   - Download raw data and clinical metadata from GEO or TCGA-GDC.
   - Perform platform-specific normalization (RMA, quantile normalization, or `rlog`) and differential expression analysis (`limma` or `DESeq2`).
   - Save high-resolution vector visualizations (`.svg`), including Boxplot, PCA, UMAP, Mean-Variance, Mean-Difference, Quantile-Quantile plot, Volcano plot, Adjusted p-value histogram, Venn diagram, Heatmap, Density plot, and Scatter plot.
   - Export a Z-score normalized `*_ML_Ready.csv` file inside each dataset folder *(for `GSE61710`, peritumoral samples are excluded prior to Z-score scaling to maintain a strict GBM vs. Normal binary classification)*.
   - *(TCGA only)* Run Kaplan-Meier survival analysis using both MaxStat optimal cutpoints and a pre-specified median split, compute the multivariate **GBM miRNA Prognostic Index (GMPI)**, quantify discrimination via **Harrell's C-index** (with 1,000 bootstrap resamples) and **1-, 2-, and 3-year time-dependent ROC-AUCs** (`timeROC`), fit an expanded multivariate Cox model adjusted for clinical covariates (`age`, `sex`, and `radiation`), and benchmark clinical-only versus combined models via **Likelihood-Ratio Tests (LRT)** and **1-year Decision Curve Analysis** (`dcurves`).

### Step 2: Cross-Cohort Meta-Analysis (R)
1. While still inside the `Main.Rproj` RStudio session, open and run the R script inside the `Meta Analysis/` folder.
2. This script reads the differential expression results from all seven dataset folders, standardizes miRNA nomenclature, and generates:
   - Overlap visualizations (`MetaAnalysis_All_VennDiagram` and `MetaAnalysis_All_UpSetPlot` in `.png` and `.svg`).
   - A list of universally detected miRNAs (`Conserved_All_miRNAs.csv`).
   - Fisher's Combined Probability Test results (`MetaAnalysis_Fisher_Results.csv`) and a labeled meta-volcano plot (`MetaAnalysis_Custom_Volcano_Labeled`).

### Step 3: Unified Machine Learning, LODO Validation & Explainable AI (Python)
Run `ML_pipeline.py` from the root directory:
```bash
python "ML_pipeline.py"
```
This unified, fully deterministic script executes four sequential modules:
1. **Automated Cohort Merging**: Recursively scans all dataset subfolders for `*_ML_Ready.csv` files, tags each sample with its source cohort (`Dataset` column), performs an inner join to retain only universally detected (conserved) miRNAs across all seven cohorts, and outputs `MASTER_GBM_DATASET.csv`.
2. **Internal 10×5-Fold CV & Ablation Study**: Runs 10 iterations of balanced random undersampling with nested 5-fold stratified cross-validation (50 total folds, fixed random seeds) comparing five feature selection configurations (`RF_Alone`, `LASSO_Alone`, `SVM_Alone`, `Consensus_Unweighted`, and `Consensus_Weighted` using network degree centrality). Reports **Mean ± SD** and **95% confidence intervals** across folds, performs fold-level paired significance testing (**Nadeau–Bengio corrected repeated-CV *t*-test** and **paired Wilcoxon signed-rank test**) against `Consensus_Weighted`, and exports `Ablation_Study_Results.csv`.
3. **7-Fold Leave-One-Dataset-Out (LODO) & Pooled LODO Validation**: Iteratively trains on 6 cohorts (with nested topologically-weighted feature selection and 10-iteration balanced undersampling) and evaluates on the 7th unseen cohort using a fixed 0.50 probability threshold. Computes per-cohort AUCs, cohort sample-size weights (`Weight_%`), unweighted and sample-size-weighted mean AUCs, and a **Pooled LODO AUC** with **95% bootstrap confidence intervals** (1,000 resamples). Exports `LODO_Validation_Results.csv` and `ML_Outputs/LODO_ROC_Curves` (`.png` and `.svg`).
4. **Explainable AI (SHAP) & High-Definition ROC**: Saves `SHAP_Summary_Plot`, individual `SHAP_Dependence_*` plots for the 4-miRNA panel, and `High_Def_ROC_Curve` (600 DPI `.png` and vector `.svg`) inside `ML_Outputs/`.

---

## Notes on Network Access
Automated downloads from NCBI GEO and the NCI GDC portal (`GEOquery` and `TCGAbiolinks`) require a stable internet connection. If you encounter connection timeouts or regional network restrictions, you can route R through a local proxy or manually download the raw archives into the corresponding dataset folder before running the scripts.

## License
This project is licensed under the terms specified in the `LICENSE` file.
