# ==============================================================================
# 1. SETUP AND PACKAGE INSTALLATION
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages_needed <- c(
  "dplyr", "TCGAbiolinks", "edgeR", "limma", "DESeq2", 
  "ggplot2", "EnhancedVolcano", "pheatmap", "gplots", 
  "miRBaseConverter", "umap", "survival", "survminer", 
  "readr", "ggpubr", "maxstat","ggrepel"
)

for (p in packages_needed) {
  if (!require(p, character.only = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}
for (p in packages_needed) {
  library(p, character.only = TRUE)
}

# ==============================================================================
# 2. DIRECTORY SETUP & DATA ACQUISITION (TCGAbiolinks)
# ==============================================================================
# Anchored by Main.Rproj at the root folder.
base_dir <- "TCGA_RNAseq"
meta_dir <- "TCGA_RNAseq/Meta"

if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)

message("Querying TCGA-GBM miRNA data...")
query_TCGA <- GDCquery(project = "TCGA-GBM",
                       data.category = "Transcriptome Profiling",
                       data.type = "miRNA Expression Quantification",
                       workflow.type = "BCGSC miRNA Profiling",
                       access = "open")

query_result <- getResults(query_TCGA)
write.table(query_result, file = file.path(base_dir, "all_samples_of_GBM.txt"), quote = FALSE, sep = "\t")

# Isolate Normal and Tumor samples
normal_sample_list <- subset(query_result, sample_type == "Solid Tissue Normal" & experimental_strategy == "miRNA-Seq")
write.table(normal_sample_list, file = file.path(base_dir, "all_normal_samples_of_GBM.txt"), quote = FALSE, sep = "\t")

Tumor_sample_list <- subset(query_result, sample_type == "Primary Tumor" & experimental_strategy == "miRNA-Seq")
write.table(Tumor_sample_list, file = file.path(base_dir, "all_tumor_samples_of_GBM.txt"), quote = FALSE, sep = "\t")

# Combine the sample barcodes
all_GBM_samples <- c(normal_sample_list$cases, Tumor_sample_list$cases)
write.table(all_GBM_samples, file = file.path(base_dir, "all_GBM_access_status.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

selected_sample_vec <- read.table(file.path(base_dir, "all_GBM_access_status.txt"))[[1]]
query_TCGA_2 <- GDCquery(project = "TCGA-GBM",
                         data.category = "Transcriptome Profiling",
                         data.type = "miRNA Expression Quantification",
                         workflow.type = "BCGSC miRNA Profiling",
                         barcode = selected_sample_vec)

GDCdownload(query_TCGA_2)

# Prepare the data by reading from the new Samples folder
data_prepare <- GDCprepare(query_TCGA_2, 
                           save = TRUE, 
                           save.filename = file.path(base_dir, "data.rda"))

write.table(data_prepare, file.path(base_dir, "Raw_data.txt"), quote=FALSE, sep="\t", row.names = FALSE)

miRNA_col <- data_prepare[, 1, drop = FALSE] 
count_cols <- data_prepare[, grepl("^read_count_", colnames(data_prepare))]
df_counts <- cbind(miRNA_col, count_cols)
write.table(df_counts, file.path(base_dir, "GBM_miRNA_counts_only.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

miRNA_names <- as.character(df_counts[[1]])
rownames(df_counts) <- miRNA_names
df_counts <- df_counts[,-1]
write.table(df_counts, file.path(base_dir, "df_counts.txt"), quote = FALSE, sep = "\t")

# ==============================================================================
# 3. PREPROCESSING & DIFFERENTIAL EXPRESSION (DESeq2)
# ==============================================================================
counts <- read.table(file.path(base_dir, "GBM_miRNA_counts_only.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(counts) <- sub("^read_count_", "", colnames(counts))
rownames(counts) <- counts[,1]
counts <- counts[,-1]

# Create group vector based on TCGA barcode sample type (11 = Solid Tissue Normal, 01 = Primary Solid Tumor)
sample_types <- substr(colnames(counts), 14, 15)
group <- factor(ifelse(sample_types == "11", "Normal", ifelse(sample_types == "01", "Tumor", "Other")))

y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

write.table(cbind(miRNA_ID = rownames(y), y$counts), file = file.path(base_dir, "GBM_miRNA_filtered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

data_filtered <- read.table(file.path(base_dir, "GBM_miRNA_filtered.txt"), header = TRUE, row.names = 1, check.names = FALSE)
selected_sample <- colnames(data_filtered)
sample_type <- ifelse(substr(selected_sample, 14, 15) == "11", "Normal", ifelse(substr(selected_sample, 14, 15) == "01", "Tumor", "Other"))

sample_info <- data.frame(selected_sample = selected_sample, sample_type = sample_type, stringsAsFactors = FALSE)
write.table(sample_info, file.path(base_dir, "MetaData.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

data_filtered2 <- as.data.frame(lapply(data_filtered, as.numeric), row.names = rownames(data_filtered))

# DESeq2 Execution
d1 <- DESeqDataSetFromMatrix(countData = data_filtered2, colData = sample_info, design = ~sample_type)
data_rlog <- rlogTransformation(d1)
data_rlog_df <- data.frame(assay(data_rlog))
write.table(data_rlog_df, file.path(base_dir, "rLog_Data.txt"), quote = FALSE, sep = "\t")

d2 <- estimateSizeFactors(d1)
d3 <- DESeq(d2)

deseq_result <- results(d3, pAdjustMethod = "fdr")
final_result <- data.frame(deseq_result)
write.table(final_result, file.path(base_dir, "DEG.txt"), quote = FALSE, sep = "\t")

final_result_sort <- final_result[order(final_result$padj, decreasing = FALSE),]
write.table(final_result_sort, file.path(base_dir, "Final_Result_sorted.txt"), quote = FALSE, sep = "\t")

data_rlog_df2 <- as.matrix(read.table(file.path(base_dir, "rLog_Data.txt"), sep = "\t", header = TRUE, row.names = 1))

# Define unified top 10 miRNAs based on statistical significance (DESeq2 padj)
top_10_miRNAs <- rownames(final_result_sort)[1:10]
top_50_miRNAs <- rownames(final_result_sort)[1:50]
# ==============================================================================
# 4. VISUALIZATIONS (PNG, PDF & SVG Output) 
# ==============================================================================
message("Generating high-resolution plots in SVG format...")

group_code <- substr(colnames(data_rlog_df2), 14, 15)
conditions <- factor(ifelse(group_code == "11", "Normal", "Tumor"), levels = c("Normal", "Tumor"))
final_result2 <- final_result_sort
colnames(final_result2)[colnames(final_result2) == "padj"] <- "Adj.pValue"
# ------------------------------------------------------------------------------
# 1a. Boxplot (Raw vs Normalized Data)
# ------------------------------------------------------------------------------
# --- Plot 1: 50-Sample Subset ---
draw_boxplots_subset <- function() {
  par(mfrow=c(1,2), mar=c(7, 4, 4, 2) + 0.1)
  set.seed(42) 
  subset_cols <- sample(1:ncol(data_filtered2), 50)
  box_colors <- ifelse(conditions[subset_cols] == "Normal", "blue", "red")
  
  raw_exprs <- log2(data.matrix(data_filtered2) + 1)
  boxplot(raw_exprs[, subset_cols], main="Raw Data (Subset)", outline=FALSE, las=2, col=box_colors, xaxt='n', ylab="Log2(Counts + 1)")
  legend("topright", legend=c("Normal", "GBM"), fill=c("blue", "red"))
  
  norm_exprs <- data.matrix(data_rlog_df2)
  boxplot(norm_exprs[, subset_cols], main="Normalized Data (Subset)", outline=FALSE, las=2, col=box_colors, xaxt='n', ylab="rlog Expression")
  legend("topright", legend=c("Normal", "GBM"), fill=c("blue", "red"))
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

svg(file.path(base_dir, "Plot_01a_Boxplot_Subset.svg"), width=12, height=6); draw_boxplots_subset(); dev.off()

# --- Plot 2: All Samples ---
draw_boxplots_all <- function() {
  par(mfrow=c(1,2), mar=c(7, 4, 4, 2) + 0.1)
  box_colors <- ifelse(conditions == "Normal", "blue", "red")
  
  raw_exprs <- log2(data.matrix(data_filtered2) + 1)
  boxplot(raw_exprs, main="Raw Data (All Samples)", outline=FALSE, las=2, col=box_colors, xaxt='n', ylab="Log2(Counts + 1)")
  legend("topright", legend=c("Normal", "GBM"), fill=c("blue", "red"))
  
  norm_exprs <- data.matrix(data_rlog_df2)
  boxplot(norm_exprs, main="Normalized Data (All Samples)", outline=FALSE, las=2, col=box_colors, xaxt='n', ylab="rlog Expression")
  legend("topright", legend=c("Normal", "GBM"), fill=c("blue", "red"))
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

svg(file.path(base_dir, "Plot_01a_Boxplot_All.svg"), width=12, height=6); draw_boxplots_all(); dev.off()
# ------------------------------------------------------------------------------
# 1b. Sample Library Size Barplots (Quality Control)
# ------------------------------------------------------------------------------
draw_chunked_barplots <- function() {
  sum_df <- data.frame(name_list = colnames(data_rlog_df2), colsums_list = colSums(data_rlog_df2))
  chunk_size <- 20
  split_data <- split(sum_df, (seq(nrow(sum_df)) - 1) %/% chunk_size)
  
  for (i in seq_along(split_data)) {
    p <- ggplot(split_data[[i]], aes(x = name_list, y = colsums_list)) +
      geom_bar(stat = "identity", fill = "tomato") +
      labs(title = paste("QC Barplot (Total rlog Expression) - Samples", 
                         min((i - 1) * chunk_size + 1, nrow(sum_df)), "to", 
                         min(i * chunk_size, nrow(sum_df))),
           x = "Samples", y = "Total Expression") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
    
    file_base <- file.path(base_dir, paste0("Plot_01b_QC_Barplot_Chunk_", i))
    ggsave(paste0(file_base, ".svg"), plot = p, width = 10, height = 6, units="in", dpi=300)
  }
}
draw_chunked_barplots()

# ------------------------------------------------------------------------------
# 2. PCA Plot
# ------------------------------------------------------------------------------
pca <- prcomp(t(data_rlog_df2), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = conditions)

pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition)) +
  geom_point(size=3) + theme_bw() + ggtitle("PCA: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red"), labels=c("Normal"="Normal", "Tumor"="GBM"))

ggsave(file.path(base_dir, "Plot_02_PCA.svg"), plot=pca_plot, width=7, height=6, units="in", dpi=300)

# ------------------------------------------------------------------------------
# 3. UMAP Plot
# ------------------------------------------------------------------------------
custom_umap <- umap.defaults; custom_umap$n_neighbors <- 15 
umap_out <- umap(t(data_rlog_df2), config=custom_umap)
umap_df <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Condition = conditions)

umap_plot <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=Condition)) +
  geom_point(size=3, alpha=0.6) + theme_bw() + ggtitle("UMAP: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red"), labels=c("Normal"="Normal", "Tumor"="GBM"))

ggsave(file.path(base_dir, "Plot_03_UMAP.svg"), plot=umap_plot, width=7, height=6, units="in", dpi=300)

# ------------------------------------------------------------------------------
# 4. Mean-Variance Trend (limma-voom approach)
# ------------------------------------------------------------------------------
draw_mv <- function() {
  y_norm <- calcNormFactors(y)
  design <- model.matrix(~ sample_type, data = sample_info)
  v <- voom(y_norm, design, plot = FALSE)
  fit_voom <- lmFit(v, design)
  fit_voom <- eBayes(fit_voom)
  plotSA(fit_voom, main = "Mean-Variance Trend")
}
svg(file.path(base_dir, "Plot_04_Mean_Variance_Trend.svg"), width=8, height=6); draw_mv(); dev.off()

# ------------------------------------------------------------------------------
# 5. Q-Q Plot
# ------------------------------------------------------------------------------
draw_qq <- function() {
  qqnorm(final_result$stat, main="Q-Q Plot", xlab="Theoretical Quantiles", ylab="Sample Quantiles", pch=16, col="darkgray", cex=0.6)
  qqline(final_result$stat, col="red", lwd=2)
}
svg(file.path(base_dir, "Plot_05_QQ_Plot.svg"), width=8, height=6); draw_qq(); dev.off()

# ------------------------------------------------------------------------------
# 6. Mean-Difference (MD / MA) Plots
# ------------------------------------------------------------------------------
draw_ma <- function() { plotMA(deseq_result, main="MA Plot (DESeq2: GBM vs Normal)", ylim=c(-5,5)) }
svg(file.path(base_dir, "Plot_06_MA_Plot_DESeq2.svg"), width=8, height=6); draw_ma(); dev.off()

draw_md_voom <- function() {
  y_norm <- calcNormFactors(y)
  design <- model.matrix(~ sample_type, data = sample_info)
  v <- voom(y_norm, design, plot = FALSE)
  fit_voom <- lmFit(v, design)
  fit_voom <- eBayes(fit_voom)
  plotMD(fit_voom, column = 2, status = decideTests(fit_voom), main = "MD Plot (Voom): GBM vs Normal")
}
svg(file.path(base_dir, "Plot_06_MD_Plot_Voom.svg"), width=8, height=6); draw_md_voom(); dev.off()

# ------------------------------------------------------------------------------
# 7. Volcano Plot
# ------------------------------------------------------------------------------
volcano <- EnhancedVolcano(final_result2,
                           lab = rownames(final_result2), x = 'log2FoldChange', y = 'Adj.pValue',
                           selectLab = top_10_miRNAs,
                           title = 'RNA-seq: GBM vs Normal',
                           pCutoff = 0.05, FCcutoff = 1.0,
                           pointSize = 3.0, labSize = 3.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5)
ggsave(file.path(base_dir, "Plot_07_Volcano.svg"), plot = volcano, width = 12, height = 10, units="in", dpi=300)

# ------------------------------------------------------------------------------
# 8. Adjusted p value histogram
# ------------------------------------------------------------------------------
draw_hist <- function() { hist(final_result$padj, col = "lightblue", main = "Histogram of Adj P-values", xlab = "Adjusted P-value", breaks = 50) }
svg(file.path(base_dir, "Plot_08_Adj_PVal_Hist.svg"), width=8, height=6); draw_hist(); dev.off()

# ------------------------------------------------------------------------------
# 9. Venn Diagram
# ------------------------------------------------------------------------------
draw_venn <- function() {
  y_norm <- calcNormFactors(y)
  design <- model.matrix(~ 0 + sample_type, data = sample_info)
  colnames(design) <- gsub("sample_type", "", colnames(design))
  colnames(design)[colnames(design) == "Tumor"] <- "GBM" 
  
  v <- voom(y_norm, design, plot = FALSE)
  fit_voom <- lmFit(v, design)
  cont.matrix <- makeContrasts(GBMvsNormal = GBM - Normal, levels = design)
  fit2_voom <- contrasts.fit(fit_voom, cont.matrix)
  fit3_voom <- eBayes(fit2_voom)
  
  results_decision <- decideTests(fit3_voom, adjust.method = "fdr", p.value = 0.05)
  vennDiagram(results_decision, include = c("up", "down"), counts.col = c("red", "blue"))
}
svg(file.path(base_dir, "Plot_09_VennDiagram.svg"), width=8, height=8); draw_venn(); dev.off()

# ------------------------------------------------------------------------------
# 10. Heatmaps (Full Dataset & 5x5 Subset)
# ------------------------------------------------------------------------------
# Define the genes to plot (using top 50 based on DESeq2 significance)
top_genes_count_all <- data_rlog_df2[top_50_miRNAs , ]
normal_samples <- colnames(data_rlog_df2)[conditions == "Normal"]
tumor_samples  <- colnames(data_rlog_df2)[conditions == "Tumor"]

# ==========================================
# HEATMAP 1: ALL SAMPLES
# ==========================================
top_genes_count_full <- top_genes_count_all[, c(normal_samples, tumor_samples)]
annotation_df_full <- data.frame(
  Condition = factor(c(rep("Normal", length(normal_samples)), 
                       rep("GBM", length(tumor_samples))), 
                     levels = c("Normal", "GBM"))
)
rownames(annotation_df_full) <- colnames(top_genes_count_full)

generate_heatmap_all <- function() {
  pheatmap(top_genes_count_full, 
           scale = "row", 
           annotation_col = annotation_df_full, 
           show_colnames = FALSE,   # Keep hidden due to high sample count
           cluster_cols = FALSE, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
           main = "Top 50 Differentially Expressed miRNAs (All Samples)")
}
svg(file.path(base_dir, "Plot_10a_Heatmap_All_Samples.svg"), width = 8, height = 10); generate_heatmap_all(); dev.off()

# ==========================================
# HEATMAP 2: SUBSET (5 NORMAL VS 5 GBM)
# ==========================================
subset_normal <- normal_samples[1:5]
subset_tumor <- tumor_samples[1:5]

top_genes_count_subset <- top_genes_count_all[, c(subset_normal, subset_tumor)]

# Rename columns to clean Normal-1, Tumor-1 formats for readability
subset_names <- c(paste0("Normal-", 1:5), paste0("GBM-", 1:5))
colnames(top_genes_count_subset) <- subset_names

annotation_df_subset <- data.frame(
  Condition = factor(c(rep("Normal", 5), rep("GBM", 5)), levels = c("Normal", "GBM"))
)
rownames(annotation_df_subset) <- colnames(top_genes_count_subset)

generate_heatmap_subset <- function() {
  pheatmap(top_genes_count_subset, 
           scale = "row", 
           annotation_col = annotation_df_subset, 
           show_colnames = TRUE,   # Visible column names for the 10 samples
           cluster_cols = FALSE, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
           main = "Top 50 Differentially Expressed miRNAs (5 Normal vs 5 GBM)")
}

svg(file.path(base_dir, "Plot_10b_Heatmap_Subset.svg"), width = 8, height = 10); generate_heatmap_subset(); dev.off()

# ------------------------------------------------------------------------------
# 11. Density Plot & Adjusted P-value Hist
# ------------------------------------------------------------------------------

draw_density <- function() {
  par(mfrow=c(1,2))
  raw_exprs <- log2(data.matrix(data_filtered2) + 1)
  plotDensities(raw_exprs, group = conditions, main = "Raw Density", legend = "topright", col = c("blue", "red"))
  plotDensities(data.matrix(data_rlog_df2), group = conditions, main = "Normalized Density", legend = "topright", col = c("blue", "red"))
  par(mfrow=c(1,1))
}
svg(file.path(base_dir, "Plot_11_Density.svg"), width=12, height=6); draw_density(); dev.off()

# ------------------------------------------------------------------------------
# 12. Scatter Plot (Mean Normal vs Mean Tumor)
# ------------------------------------------------------------------------------
draw_scatter <- function() {
  normal_idx <- which(conditions == "Normal"); tumor_idx  <- which(conditions == "Tumor")
  mean_normal <- rowMeans(data_rlog_df2[, normal_idx]); mean_tumor  <- rowMeans(data_rlog_df2[, tumor_idx])
  common_genes <- intersect(rownames(final_result), rownames(data_rlog_df2))
  res_subset <- final_result[common_genes, ]; mean_n_sub <- mean_normal[common_genes]; mean_t_sub <- mean_tumor[common_genes]
  
  is_sig <- res_subset$padj < 0.05 & abs(res_subset$log2FoldChange) > 1
  is_sig[is.na(is_sig)] <- FALSE
  
  # 1. Build a dataframe for ggplot
  scatter_df <- data.frame(
    miRNA = names(mean_n_sub),
    Mean_Normal = mean_n_sub,
    Mean_GBM = mean_t_sub,
    Significance = ifelse(is_sig, "Significant", "Not Significant")
  )
  
  # 2. Create a Label column that only contains text for your top 10 miRNAs
  scatter_df$Label <- ifelse(scatter_df$miRNA %in% top_10_miRNAs, scatter_df$miRNA, "")
  
  # 3. Build the plot with ggrepel
  p <- ggplot(scatter_df, aes(x = Mean_Normal, y = Mean_GBM)) +
    geom_point(aes(color = Significance), size = 1.5, alpha = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 1) +
    geom_text_repel(aes(label = Label),
                    size = 3.5,
                    fontface = "bold",
                    box.padding = 0.8,
                    point.padding = 0.3,
                    segment.color = "black",
                    segment.size = 0.5,
                    min.segment.length = 0, # Forces the connector line to always draw
                    max.overlaps = Inf) +
    theme_bw() +
    labs(x = "Average Normal Expression (Log2)",
         y = "Average GBM Expression (Log2)",
         title = "Scatter Plot: Mean Normal vs Mean GBM") +
    theme(legend.position = "none")
  
  return(p)
}
# Generate and save the ggplot
scatter_plot <- draw_scatter()
ggsave(file.path(base_dir, "Plot_12_Scatter_Means.svg"), plot = scatter_plot, width = 8, height = 8)

# ==============================================================================
# 5. MACHINE LEARNING DATA PREPARATION
# ==============================================================================
rna_matrix <- data_rlog_df2
precursors <- rownames(rna_matrix)
mapping <- miRNA_PrecursorToMature(precursors)

final_rows <- c()
final_names <- c()

for (i in 1:nrow(mapping)) {
  orig_name <- mapping$OriginalName[i]
  mat1 <- mapping$Mature1[i]
  mat2 <- mapping$Mature2[i]
  
  if (orig_name %in% rownames(rna_matrix)) {
    expr_vals <- rna_matrix[orig_name, ]
    if (!is.na(mat1) && mat1 != "") {
      final_rows <- rbind(final_rows, expr_vals)
      final_names <- c(final_names, mat1)
    }
    if (!is.na(mat2) && mat2 != "") {
      final_rows <- rbind(final_rows, expr_vals)
      final_names <- c(final_names, mat2)
    }
  }
}

rownames(final_rows) <- final_names
colnames(final_rows) <- colnames(rna_matrix)
ml_matrix <- t(final_rows)
z_score_matrix <- scale(ml_matrix)
ml_df <- as.data.frame(z_score_matrix, check.names = FALSE)

sample_codes <- substr(rownames(ml_df), 14, 15)
ml_df$Target <- ifelse(sample_codes == "11", 0, 1)
ml_df$Sample_ID <- rownames(ml_df)
mirna_cols <- setdiff(colnames(ml_df), c("Sample_ID", "Target"))
ml_df <- ml_df[, c("Sample_ID", mirna_cols, "Target")]
ml_df <- ml_df[, !duplicated(colnames(ml_df))]

write.csv(ml_df, file.path(base_dir, "TCGA_RNASeq_ML_Ready.csv"), row.names = FALSE)

# ==============================================================================
# 6. KAPLAN-MEIER SURVIVAL ANALYSIS & COX PROPORTIONAL HAZARDS (NOVELTY)
# ==============================================================================
message("Downloading clinical data directly from TCGA...")

# Automatically fetch the clinical data from the GDC portal
clinical_data <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")

# Save a local copy in your Meta directory for record-keeping
write.table(clinical_data, file.path(meta_dir, "TCGA-GBM_clinical_auto.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

expr_path <- file.path(base_dir, "TCGA_RNASeq_ML_Ready.csv")
expr_data <- read.csv(expr_path, stringsAsFactors = FALSE)

# Map the TCGAbiolinks column names to match your survival pipeline
clinical_clean <- clinical_data %>%
  select(case_submitter_id = submitter_id,           
         vital_status = vital_status,          
         days_to_death = days_to_death,        
         days_to_last_follow_up = days_to_last_follow_up,
         age = age_at_index) %>%
  mutate(across(everything(), ~na_if(as.character(.), "'--"))) %>%
  mutate(across(everything(), ~na_if(as.character(.), "not reported"))) %>%
  unique()

clinical_clean$days_to_death <- as.numeric(clinical_clean$days_to_death)
clinical_clean$days_to_last_follow_up <- as.numeric(clinical_clean$days_to_last_follow_up)
clinical_clean$age <- as.numeric(clinical_clean$age)

clinical_clean <- clinical_clean %>%
  mutate(OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
         OS_status = ifelse(vital_status == "Dead", 1, 0)) %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(age))

expr_tumor <- expr_data %>% filter(Target == 1)
expr_tumor$case_submitter_id <- gsub("\\.", "-", substr(expr_tumor$Sample_ID, 1, 12))
expr_tumor_dedup <- expr_tumor %>% distinct(case_submitter_id, .keep_all = TRUE)
clinical_clean_dedup <- clinical_clean %>% distinct(case_submitter_id, .keep_all = TRUE)
merged_data <- inner_join(expr_tumor_dedup, clinical_clean_dedup, by = "case_submitter_id")

target_mirna <- c("hsa.miR.769.3p", "hsa.miR.485.3p", "hsa.miR.137", "hsa.miR.765")

for (mirna in target_mirna) {
  if (mirna %in% colnames(merged_data)) {
    res.cut <- surv_cutpoint(merged_data, time = "OS_time", event = "OS_status", variables = mirna, minprop = 0.2)
    res.cat <- surv_categorize(res.cut)
    merged_data$expr_group <- factor(res.cat[[mirna]], levels = c("low", "high"))
    
    surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OS_status)
    fit <- survfit(surv_object ~ expr_group, data = merged_data)
    
    p <- ggsurvplot(fit, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("#2E9FDF", "#E7B800"), title = paste("Survival (Optimal Cutpoint):", mirna))
    
    base_name <- file.path(base_dir, paste0("Plot_13_KM_Optimal_Cutpoint_", mirna))
    svg(filename = paste0(base_name, ".svg"), width = 8, height = 7); print(p); dev.off()
  }
}

# Exo-GBM Prognostic Index (EGPI) Custom Risk Score
formula_str <- paste0("Surv(time = OS_time, event = OS_status) ~ `", paste(target_mirna, collapse = "` + `"), "` + age")
cox_formula <- as.formula(formula_str)
cox_model <- coxph(cox_formula, data = merged_data)

merged_data$EGPI <- predict(cox_model, type = "risk")
res.cut_egpi <- surv_cutpoint(merged_data, time = "OS_time", event = "OS_status", variables = "EGPI", minprop = 0.2)
res.cat_egpi <- surv_categorize(res.cut_egpi)
merged_data$EGPI_Group <- factor(res.cat_egpi$EGPI, levels = c("low", "high"))

fit_opt_egpi <- survfit(Surv(time = OS_time, event = OS_status) ~ EGPI_Group, data = merged_data)
p_opt_egpi <- ggsurvplot(fit_opt_egpi, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE, palette = c("#2E9FDF", "#E7B800"), title = "Survival: Exo-GBM Prognostic Index (EGPI)")
svg(filename = file.path(base_dir, "Plot_14_KM_Multivariate_EGPI_Optimal_Cutpoint.svg"), width = 8, height = 7); print(p_opt_egpi); dev.off()