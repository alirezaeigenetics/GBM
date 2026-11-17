# ==============================================================================
# SETUP AND PACKAGE INSTALLATION
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages_needed <- c("limma", "umap", "EnhancedVolcano", "pheatmap","ggrepel", "ggplot2", "statmod", "GEOquery")
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
# 1. DIRECTORY SETUP & AUTOMATED DATA DOWNLOAD
# ==============================================================================
# Anchored by Main.Rproj at the root folder.
gpr_path <- "GSE93850/Samples"

if (!dir.exists(gpr_path)) {
  dir.create(gpr_path, recursive = TRUE)
}

# Download and extract the raw TXT/GPR files if the folder is empty
if (length(list.files(gpr_path, pattern = "\\.txt(\\.gz)?|\\.gpr(\\.gz)?$", ignore.case = TRUE)) == 0) {
  message("Downloading raw data files from GEO (GSE93850)... This may take a moment.")
  
  tryCatch({
    getGEOSuppFiles("GSE93850", makeDirectory = FALSE, baseDir = gpr_path)
    tar_file <- file.path(gpr_path, "GSE93850_RAW.tar")
    
    if (file.exists(tar_file)) {
      untar(tar_file, exdir = gpr_path)
      file.remove(tar_file) 
    }
  }, error = function(e) {
    stop("\n[Network Error]: Automated download from NCBI GEO failed (Possible IP block or timeout).\n",
         "Please download the dataset manually:\n",
         "1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93850\n",
         "2. Download 'GSE93850_RAW.tar' (bottom of the page)\n",
         "3. Extract the files into the 'GSE93850/Samples' folder.\n",
         "4. Re-run this script.")
  })
}

# ==============================================================================
# 2. DYNAMIC METADATA MAPPING & NORMALIZATION
# ==============================================================================
message("Fetching clinical metadata from GEO to map conditions...")
gse <- getGEO("GSE93850", GSEMatrix = TRUE,destdir = "GSE93850")
pdata <- pData(gse[[1]])

# Target the 'title' column which safely contains the "Normal" / "Tumor" naming convention
pdata$Condition <- ifelse(grepl("normal", pdata$title, ignore.case = TRUE), "Normal", "Tumor")

extracted_files <- list.files(gpr_path, pattern = "\\.txt(\\.gz)?|\\.gpr(\\.gz)?$", ignore.case = TRUE)

# 1. Try to map from GEO metadata based on GSM ID
geo_conditions <- pdata$Condition[match(substr(extracted_files, 1, 10), rownames(pdata))]

# 2. Fallback for custom file names if GEO match fails
final_conditions <- ifelse(!is.na(geo_conditions), geo_conditions, 
                           ifelse(grepl("Normal", extracted_files, ignore.case = TRUE), "Normal", 
                                  ifelse(grepl("Tumor|GBM", extracted_files, ignore.case = TRUE), "Tumor", NA)))

targets <- data.frame(
  FileName = extracted_files,
  Condition = final_conditions
)

missing_files <- targets$FileName[!file.exists(file.path(gpr_path, targets$FileName))]
if(length(missing_files) > 0) stop("Missing files:\n", paste(missing_files, collapse="\n"))
message(paste("All", nrow(targets), "files successfully mapped. Proceeding to cleaning..."))

# Function to fix malformed Agilent GPR file headers
clean_gpr_header <- function(filename, folder_path) {
  full_path <- file.path(folder_path, filename)
  lines <- readLines(full_path, warn = FALSE, encoding = "latin1")
  bad_line_index <- grep("^ImageFiles=", lines)
  if (length(bad_line_index) > 0) {
    lines[bad_line_index] <- "ImageFiles=CleanedByR.tif" 
    writeLines(lines, full_path)
    return(TRUE)
  }
  return(FALSE)
}

# Apply the header cleaner to all target files before reading
invisible(sapply(targets$FileName, clean_gpr_header, folder_path = gpr_path))

# Read and normalize Agilent data
raw_data <- read.maimages(files = targets$FileName, path = gpr_path, source = "generic", columns = list(E = "F547 Median", Eb = "B547 Median"), annotation = c("Block", "Column", "Row", "Name", "ID"), green.only = TRUE)
colnames(raw_data) <- paste(targets$Condition, ave(as.character(targets$Condition), targets$Condition, FUN = seq_along), sep="_")

bg_corrected <- backgroundCorrect(raw_data, method="normexp", offset=50)
norm_data <- normalizeBetweenArrays(bg_corrected, method="quantile")

if("Name" %in% colnames(norm_data$genes)){
  valid_probes <- norm_data$genes$Name != "" & !is.na(norm_data$genes$Name)
  norm_data_filt <- norm_data[valid_probes, ]
  norm_data_avg <- avereps(norm_data_filt, ID=norm_data_filt$genes$Name)
} else {
  norm_data_avg <- avereps(norm_data, ID=norm_data$genes$ID)
}

expr_matrix_final <- norm_data_avg$E

# ==============================================================================
# 3. FILTERING & DIFFERENTIAL EXPRESSION (limma)
# ==============================================================================
# 1. Forgiving human miRNA filter
human_indices <- grep("hsa", rownames(expr_matrix_final), ignore.case = TRUE)
expr_matrix_final <- expr_matrix_final[human_indices, ]

if(nrow(expr_matrix_final) == 0) stop("Error: No human miRNAs (hsa) found. Check annotation mapping.")

# 2. Identify and keep ONLY definitively mapped samples
valid_idx <- which(!is.na(targets$Condition) & targets$Condition %in% c("Normal", "Tumor"))

targets <- targets[valid_idx, , drop = FALSE]
expr_matrix_final <- expr_matrix_final[, valid_idx, drop = FALSE]

if(ncol(expr_matrix_final) == 0) stop("Error: No samples mapped to valid conditions.")

# 3. Build the exact-matching design matrix
condition_factor <- factor(targets$Condition, levels=c("Normal", "Tumor"))
design <- model.matrix(~0 + condition_factor)
colnames(design) <- c("Normal", "Tumor")

# 4. Fit the limma model
fit <- lmFit(expr_matrix_final, design)
cont.matrix <- makeContrasts(TumorVsNormal = Tumor - Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2, trend=TRUE)

# 5. Extract results
results <- topTable(fit3, adjust="fdr", number=Inf, sort.by="B")
write.csv(results, "GSE93850/GSE93850_Tumor_vs_Normal.csv", row.names = TRUE)

sig_degs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
write.csv(sig_degs, "GSE93850/GSE93850_Significant_DEGs.csv", row.names = TRUE)

# Define unified top 10 miRNAs based on statistical significance
top_10_miRNAs <- rownames(results)[1:10]

# ==============================================================================
# 4. VISUALIZATIONS 
# ==============================================================================
draw_boxplots <- function() {
  par(mfrow=c(1,2), mar=c(10, 4, 4, 2) + 0.1)
  
  color_map <- c("Normal"="blue", "Tumor"="red")
  sample_colors <- color_map[as.character(targets$Condition)]
  
  boxplot(log2(raw_data$E[, valid_idx]), main="Raw Data", las=2, outline=FALSE, col=sample_colors)
  boxplot(expr_matrix_final, main="Normalized Data", las=2, outline=FALSE, col=sample_colors)
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

svg("GSE93850/Plot_01_QC_Boxplots.svg", width=12, height=6)
draw_boxplots()
dev.off()

pca <- prcomp(t(expr_matrix_final), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = targets$Condition)
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition)) +
  geom_point(size=5) + theme_bw() + ggtitle("PCA: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red"), labels=c("Normal"="Normal", "Tumor"="GBM"))
ggsave("GSE93850/Plot_02_PCA.svg", plot=pca_plot, width=7, height=6, units="in", dpi=300)

custom_umap <- umap::umap.defaults
custom_umap$n_neighbors <- min(10, ncol(expr_matrix_final) - 1)
umap_out <- umap::umap(t(expr_matrix_final), config=custom_umap)
umap_df <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Condition = targets$Condition)
umap_plot <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=Condition)) +
  geom_point(size=4) + theme_bw() + ggtitle("UMAP: GBM vs Normal") +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red"), labels=c("Normal"="Normal", "Tumor"="GBM"))
ggsave("GSE93850/Plot_03_UMAP.svg", plot=umap_plot, width=7, height=6, units="in", dpi=300)

svg("GSE93850/Plot_04_Mean_Variance.svg", width=8, height=6)
plotSA(fit3, main="Mean-variance Trend", xlab="Average log-expression", ylab="sqrt(sigma)")
dev.off()

svg("GSE93850/Plot_05_QQ_Plot.svg", width=8, height=6)
qqt(fit3$t, df = fit3$df.total + fit3$df.prior, main="Q-Q Plot", xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(0,1, col="red")
dev.off()

svg("GSE93850/Plot_06_MD_Plot.svg", width=8, height=6)
plotMD(fit3, column = 1, status = decideTests(fit3), main = "MD Plot: GBM vs Normal", xlab = "Average log-expression", ylab = "log-fold-change", hl.col = c("blue", "red"))
dev.off()

volcano <- EnhancedVolcano(results, 
                           lab = rownames(results), 
                           x = 'logFC', 
                           y = 'adj.P.Val', 
                           selectLab = top_10_miRNAs,
                           title = 'GSE93850: GBM vs Normal', 
                           pCutoff = 0.05, 
                           FCcutoff = 1.0, 
                           pointSize = 3.0, 
                           labSize = 3.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5)
ggsave("GSE93850/Plot_07_Volcano.svg", plot=volcano, width=12, height=10, units="in", dpi=300)

svg("GSE93850/Plot_08_Adj_PVal_Hist.svg", width=8, height=6)
hist(results$adj.P.Val, col = "lightblue", main = "Histogram of Adj P-values", xlab = "Adjusted P-value", ylab = "Frequency", breaks = 20) 
dev.off()

results_decision <- decideTests(fit3, adjust.method="fdr", p.value=0.05)
svg("GSE93850/Plot_09_VennDiagram.svg", width=8, height=8)
colnames(results_decision) <- gsub("Tumor", "GBM", colnames(results_decision))
vennDiagram(results_decision, include=c("up", "down"), counts.col=c("red", "blue"))
dev.off()

top_50 <- rownames(results)[1:min(50, nrow(results))]
heatmap_data <- expr_matrix_final[top_50, ]
svg("GSE93850/Plot_10_Heatmap.svg", width=8, height=10)
pheatmap(heatmap_data, scale = "row", main = "Top 50 Differentially Expressed miRNAs", annotation_col = data.frame(Condition = targets$Condition, row.names=colnames(heatmap_data)), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

svg("GSE93850/Plot_11_Density_Plot.svg", width=12, height=6)
par(mfrow=c(1,2))
color_map <- c("Normal"="blue", "Tumor"="red")
density_colors <- color_map[unique(as.character(targets$Condition))]

plotDensities(log2(raw_data$E[, valid_idx]), group = targets$Condition, main = "Raw Density", legend = "topright", col = density_colors) 
plotDensities(expr_matrix_final, group = targets$Condition, main = "Normalized Density", legend = "topright", col = density_colors) 
par(mfrow=c(1,1))
dev.off()

# Means calculation safely handling Single-Sample Condition
normal_idx <- which(targets$Condition == "Normal")
tumor_idx  <- which(targets$Condition == "Tumor")
mean_normal <- if(length(normal_idx) > 1) rowMeans(expr_matrix_final[, normal_idx]) else expr_matrix_final[, normal_idx]
mean_tumor <- if(length(tumor_idx) > 1) rowMeans(expr_matrix_final[, tumor_idx]) else expr_matrix_final[, tumor_idx]

sig_genes <- rownames(results)[results$adj.P.Val < 0.05 & abs(results$logFC) > 1]
is_sig <- rownames(expr_matrix_final) %in% sig_genes

svg("GSE93850/Plot_12_Scatter_Means.svg", width=8, height=8)
plot(mean_normal, mean_tumor, pch = 16, cex = 0.5, col = "grey", xlab = "Average Normal Expression (Log2)", ylab = "Average GBM Expression (Log2)", main = "Scatter Plot: Mean Normal vs Mean GBM")
abline(0, 1, col = "blue", lwd = 2)
points(mean_normal[is_sig], mean_tumor[is_sig], col = "red", pch = 16, cex = 0.6)
top_10_idx <- match(top_10_miRNAs, rownames(expr_matrix_final))
text(mean_normal[top_10_idx], mean_tumor[top_10_idx], labels = top_10_miRNAs, pos = 4, cex = 0.8, font = 2, col = "black")
dev.off()

draw_scatter <- function() {
  scatter_df <- data.frame(
    miRNA = rownames(expr_matrix_final),
    Mean_Normal = mean_normal,
    Mean_GBM = mean_tumor,
    Significance = ifelse(is_sig, "Significant", "Not Significant")
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
ggsave("GSE93850/Plot_12_Scatter_Means_Connector.svg", plot = scatter_plot, width = 8, height = 8,units = "in")

# ==============================================================================
# 5. MACHINE LEARNING DATA PREPARATION
# ==============================================================================
ml_matrix <- t(expr_matrix_final)
z_score_matrix <- scale(ml_matrix)
ml_df <- as.data.frame(z_score_matrix)
ml_df$Target <- ifelse(targets$Condition == "Tumor", 1, 0)
ml_df$Sample_ID <- rownames(ml_df)
ml_df <- ml_df[, c("Sample_ID", colnames(z_score_matrix), "Target")]

write.csv(ml_df, "GSE93850/GSE93850_ML_Ready.csv", row.names = FALSE)