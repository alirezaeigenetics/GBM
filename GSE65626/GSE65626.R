# ==============================================================================
# SETUP AND PACKAGE INSTALLATION
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages_needed <- c("oligo", "limma", "umap", "EnhancedVolcano", "pheatmap","ggrepel", "ggplot2", "GEOquery")

for (p in packages_needed) {
  if (!require(p, character.only = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}

for (p in packages_needed) {
  library(p, character.only = TRUE)
}

if (!require("pd.mirna.4.0")) BiocManager::install("pd.mirna.4.0")
library(pd.mirna.4.0)

# ==============================================================================
# 1. DIRECTORY SETUP & AUTOMATED DATA DOWNLOAD
# ==============================================================================
# The working directory is anchored by Main.Rproj at the root folder.
cel_path <- "GSE65626/Samples"

# Create the folder if it does not exist
if (!dir.exists(cel_path)) {
  dir.create(cel_path, recursive = TRUE)
}

# Download and extract the raw CEL files if the folder is empty
if (length(list.files(cel_path, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE)) == 0) {
  message("Downloading raw CEL files from GEO (GSE65626)... This may take a moment.")
  
  tryCatch({
    getGEOSuppFiles("GSE65626", makeDirectory = FALSE, baseDir = cel_path)
    tar_file <- file.path(cel_path, "GSE65626_RAW.tar")
    
    if (file.exists(tar_file)) {
      untar(tar_file, exdir = cel_path)
      file.remove(tar_file) 
    }
  }, error = function(e) {
    stop("\n[Network Error]: Automated download from NCBI GEO failed (Possible IP block or timeout).\n",
         "Please download the dataset manually:\n",
         "1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65626\n",
         "2. Download 'GSE65626_RAW.tar' (bottom of the page)\n",
         "3. Extract the .CEL.gz files into the 'GSE65626/Samples' folder.\n",
         "4. Re-run this script.")
  })
}

# ==============================================================================
# 2. DYNAMIC METADATA MAPPING & NORMALIZATION
# ==============================================================================
message("Fetching clinical metadata from GEO to map Normal vs GBM conditions...")
gse <- getGEO("GSE65626", GSEMatrix = TRUE,destdir = "GSE65626")
gpl_idx <- grep("GPL19117", names(gse))

if(length(gpl_idx) == 0) stop("Error: GPL19117 metadata not found in the GEO download.")
pdata <- pData(gse[[gpl_idx]])
pdata$Condition <- ifelse(grepl("normal", pdata$title, ignore.case = TRUE), "Normal", "GBM")
extracted_files <- list.files(cel_path, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE)
geo_conditions <- pdata$Condition[match(substr(extracted_files, 1, 10), rownames(pdata))]
final_conditions <- ifelse(!is.na(geo_conditions), geo_conditions, 
                           ifelse(grepl("GBM", extracted_files, ignore.case = TRUE), "GBM", 
                                  ifelse(grepl("Normal", extracted_files, ignore.case = TRUE), "Normal", NA)))

targets <- data.frame(
  FileName = extracted_files,
  Condition = final_conditions
)

full_file_paths <- file.path(cel_path, targets$FileName)
missing_files <- targets$FileName[!file.exists(full_file_paths)]

if(length(missing_files) > 0) stop("Missing files:\n", paste(missing_files, collapse="\n"))
message(paste("All", nrow(targets), "files successfully mapped. Proceeding to normalization..."))

raw_data <- read.celfiles(filenames = full_file_paths)
norm_data <- rma(raw_data)
expr_matrix <- exprs(norm_data)
colnames(expr_matrix) <- paste(targets$Condition, ave(as.character(targets$Condition), targets$Condition, FUN = seq_along), sep="_")

gpl_id <- "GPL19117"
gpl <- getGEO(gpl_id, destdir = ".")
gpl_table <- Table(gpl)
annotated_df <- data.frame(ID = rownames(expr_matrix))
merged_anno <- merge(annotated_df, gpl_table, by = "ID", all.x = TRUE, sort = FALSE)
merged_anno <- merged_anno[match(rownames(expr_matrix), merged_anno$ID), ]
symbol_col <- ifelse("Transcript ID(Array Design)" %in% colnames(merged_anno), 
                     "Transcript ID(Array Design)", 
                     "miRNA_ID")

gene_symbols <- merged_anno[[symbol_col]]
final_symbols <- ifelse(is.na(gene_symbols) | gene_symbols == "", as.character(merged_anno$ID), as.character(gene_symbols))
rownames(expr_matrix) <- final_symbols
expr_matrix_final <- avereps(expr_matrix, ID = rownames(expr_matrix))

# ==============================================================================
# 3. FILTERING & DIFFERENTIAL EXPRESSION (limma)
# ==============================================================================
# 1. More forgiving human miRNA filter (removes the strict ^ anchor)
human_indices <- grep("hsa", rownames(expr_matrix_final), ignore.case = TRUE)
expr_matrix_final <- expr_matrix_final[human_indices, ]

if(nrow(expr_matrix_final) == 0) stop("Error: No human miRNAs (hsa) found. Check annotation mapping.")

# 2. Identify and keep ONLY definitively mapped Normal or GBM samples
valid_idx <- which(!is.na(targets$Condition) & targets$Condition %in% c("Normal", "GBM"))

# Permanently synchronize the targets and matrix (drop = FALSE prevents structural collapse)
targets <- targets[valid_idx, , drop = FALSE]
expr_matrix_final <- expr_matrix_final[, valid_idx, drop = FALSE]

if(ncol(expr_matrix_final) == 0) stop("Error: No samples mapped to Normal or GBM.")

# 3. Build the exact-matching design matrix
condition_factor <- factor(targets$Condition, levels = c("Normal", "GBM"))
design <- model.matrix(~0 + condition_factor)
colnames(design) <- c("Normal", "GBM")

# 4. Fit the limma model
fit <- lmFit(expr_matrix_final, design)
cont.matrix <- makeContrasts(GBMvsNormal = GBM - Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2, trend=TRUE)

# 5. Extract and save results
results <- topTable(fit3, adjust="fdr", number=Inf, sort.by="B")
write.csv(results, "GSE65626/GSE65626_GBM_vs_Normal.csv", row.names = TRUE)

sig_degs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
write.csv(sig_degs, "GSE65626/GSE65626_Significant_DEGs.csv", row.names = TRUE)

# Define unified top 10 miRNAs based on statistical significance for downstream plots
top_10_miRNAs <- rownames(results)[1:10]

# ==============================================================================
# 4. VISUALIZATIONS 
# ==============================================================================
draw_boxplots <- function() {
  par(mfrow=c(1,2), mar=c(10, 4, 4, 2) + 0.1) 
  boxplot(raw_data, main="Raw Data", las=2, col=c(rep("blue",3), rep("red",3)), target="core")
  boxplot(expr_matrix_final, main="Normalized Data", las=2, col=c(rep("blue",3), rep("red",3)))
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1) 
}

svg("GSE65626/Plot_01_QC_Boxplots.svg", width=12, height=6)
draw_boxplots()
dev.off()

pca <- prcomp(t(expr_matrix_final), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = targets$Condition)
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition)) +
  geom_point(size=5) + theme_bw() + ggtitle("PCA: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "GBM"="red"))
ggsave("GSE65626/Plot_02_PCA.svg", plot=pca_plot, width=7, height=6, units="in", dpi=300)

custom_umap <- umap::umap.defaults
custom_umap$n_neighbors <- min(4, ncol(expr_matrix_final) - 1) # Prevent n_neighbors > sample size
umap_out <- umap::umap(t(expr_matrix_final), config=custom_umap)
umap_df <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Condition = targets$Condition)
umap_plot <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=Condition)) +
  geom_point(size=5) + theme_bw() + ggtitle("UMAP: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "GBM"="red"))
ggsave("GSE65626/Plot_03_UMAP.svg", plot=umap_plot, width=7, height=6, units="in", dpi=300)

svg("GSE65626/Plot_04_Mean_Variance.svg", width=8, height=6)
plotSA(fit3, main="Mean-Variance Trend")
dev.off()

svg("GSE65626/Plot_05_QQ_Plot.svg", width=8, height=6)
qqt(fit3$t, df = fit3$df.total + fit3$df.prior, main="Q-Q Plot")
abline(0,1, col="red")
dev.off()

svg("GSE65626/Plot_06_MD_Plot.svg", width=8, height=6)
plotMD(fit3, column=1, status=decideTests(fit3), main="MD Plot: GBM vs Normal")
dev.off()

volcano <- EnhancedVolcano(results, 
                           lab = rownames(results), 
                           x = 'logFC', 
                           y = 'adj.P.Val', 
                           selectLab = top_10_miRNAs,
                           title = 'GSE65626: GBM vs Normal', 
                           pCutoff = 0.05, 
                           FCcutoff = 1.0, 
                           pointSize = 3.0, 
                           labSize = 3.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5)
ggsave("GSE65626/Plot_07_Volcano.svg", plot=volcano, width=12, height=10, units="in", dpi=300)

svg("GSE65626/Plot_08_Adj_PVal_Hist.svg", width=8, height=6)
hist(results$adj.P.Val, col="lightblue", main="Histogram of Adj P-values", xlab="Adjusted P-value")
dev.off()

results_decision <- decideTests(fit3, adjust.method="fdr", p.value=0.05)
colnames(results_decision) <- gsub("Tumor", "GBM", colnames(results_decision))
svg("GSE65626/Plot_09_VennDiagram.svg", width=8, height=8)
vennDiagram(results_decision, include=c("up", "down"), counts.col=c("red", "blue"))
dev.off()

top_50 <- rownames(results)[1:50]
heatmap_data <- expr_matrix_final[top_50, ]
svg("GSE65626/Plot_10_Heatmap.svg", width=8, height=10)
pheatmap(heatmap_data, scale = "row", main = "Top 50 Differentially Expressed miRNAs", annotation_col = data.frame(Condition = targets$Condition, row.names=colnames(heatmap_data)), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

svg("GSE65626/Plot_11_Density_Plot.svg", width=12, height=6)
par(mfrow=c(1,2))
plotDensities(log2(exprs(raw_data)), group = targets$Condition, main = "Raw Density", legend = "topright", col = c("blue", "red"))
plotDensities(expr_matrix_final, group = targets$Condition, main = "Normalized Density", legend = "topright", col = c("blue", "red"))
par(mfrow=c(1,1))
dev.off()

normal_idx <- which(targets$Condition == "Normal")
gbm_idx <- which(targets$Condition == "GBM")
mean_normal <- rowMeans(expr_matrix_final[, normal_idx])
mean_gbm <- rowMeans(expr_matrix_final[, gbm_idx])
log_fc_diff <- mean_gbm - mean_normal
is_significant <- abs(log_fc_diff) > 1 

svg("GSE65626/Plot_12_Scatter_Means.svg", width=8, height=8)
plot(mean_normal, mean_gbm, pch = 16, cex = 0.5, col = "grey", xlab = "Average Normal Expression (Log2)", ylab = "Average GBM Expression (Log2)", main = "Scatter Plot: Mean Normal vs Mean GBM")
abline(0, 1, col = "blue", lwd = 2)
points(mean_normal[is_significant], mean_gbm[is_significant], col = "red", pch = 16, cex = 0.6)
top_10_idx <- match(top_10_miRNAs, rownames(expr_matrix_final))
text(mean_normal[top_10_idx], mean_gbm[top_10_idx], labels = top_10_miRNAs, pos = 4, cex = 0.8, font = 2, col = "black")
dev.off()



draw_scatter <- function() {
  scatter_df <- data.frame(
    miRNA = rownames(expr_matrix_final),
    Mean_Normal = mean_normal,
    Mean_GBM = mean_gbm,
    Significance = ifelse(is_significant, "Significant", "Not Significant")
  )
  scatter_df$Label <- ifelse(scatter_df$miRNA %in% top_10_miRNAs, scatter_df$miRNA, "")
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
scatter_plot <- draw_scatter()
ggsave("GSE65626/Plot_12_Scatter_Means_Connector.svg", plot = scatter_plot, width = 8, height = 8,units = "in")

# ==============================================================================
# 5. MACHINE LEARNING DATA PREPARATION
# ==============================================================================
ml_matrix <- t(expr_matrix_final)
z_score_matrix <- scale(ml_matrix)
ml_df <- as.data.frame(z_score_matrix)
ml_df$Target <- ifelse(targets$Condition == "GBM", 1, 0)
ml_df$Sample_ID <- rownames(ml_df)
ml_df <- ml_df[, c("Sample_ID", colnames(z_score_matrix), "Target")]

write.csv(ml_df, "GSE65626/GSE65626_ML_Ready.csv", row.names = FALSE)