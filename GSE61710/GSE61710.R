# ==============================================================================
# SETUP AND PACKAGE INSTALLATION
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Added ggbio and affycoretools to resolve the deep dependency chain
packages_needed <- c("ggbio", "affycoretools", "limma", "umap", "EnhancedVolcano", "pheatmap", "ggplot2", "statmod", "AgiMicroRna", "GEOquery")

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
# The working directory is anchored by Main.Rproj at the root folder.
txt_path <- "GSE61710/Samples"

if (!dir.exists(txt_path)) {
  dir.create(txt_path, recursive = TRUE)
}

if (length(list.files(txt_path, pattern = "\\.txt(\\.gz)?$", ignore.case = TRUE)) == 0) {
  message("Downloading raw TXT files from GEO (GSE61710)... This may take a moment.")
  
  tryCatch({
    getGEOSuppFiles("GSE61710", makeDirectory = FALSE, baseDir = txt_path)
    tar_file <- file.path(txt_path, "GSE61710_RAW.tar")
    
    if (file.exists(tar_file)) {
      untar(tar_file, exdir = txt_path)
      file.remove(tar_file) 
    }
  }, error = function(e) {
    stop("\n[Network Error]: Automated download from NCBI GEO failed (Possible IP block or timeout).\n",
         "Please download the dataset manually:\n",
         "1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61710\n",
         "2. Download 'GSE61710_RAW.tar' (bottom of the page)\n",
         "3. Extract the .txt.gz files into the 'GSE61710/Samples' folder.\n",
         "4. Re-run this script.")
  })
}

# ==============================================================================
# 2. DYNAMIC METADATA MAPPING & NORMALIZATION
# ==============================================================================
message("Fetching clinical metadata from GEO to map conditions...")
gse <- getGEO("GSE61710", GSEMatrix = TRUE,destdir = "GSE61710")
pdata <- pData(gse[[1]])

pdata$Condition <- ifelse(grepl("normal", pdata$title, ignore.case = TRUE), "Normal", 
                          ifelse(grepl("peri", pdata$title, ignore.case = TRUE), "Peritumoral", "Tumor"))

extracted_files <- list.files(txt_path, pattern = "\\.txt(\\.gz)?$", ignore.case = TRUE)

# 1. Try to map from GEO metadata based on GSM ID
geo_conditions <- pdata$Condition[match(substr(extracted_files, 1, 10), rownames(pdata))]

# 2. Fallback for custom file names if GEO match fails
final_conditions <- ifelse(!is.na(geo_conditions), geo_conditions, 
                           ifelse(grepl("Normal", extracted_files, ignore.case = TRUE), "Normal", 
                                  ifelse(grepl("Peri", extracted_files, ignore.case = TRUE), "Peritumoral", 
                                         ifelse(grepl("Tumor", extracted_files, ignore.case = TRUE), "Tumor", NA))))

targets <- data.frame(
  FileName = extracted_files,
  Condition = final_conditions
)

missing_files <- targets$FileName[!file.exists(file.path(txt_path, targets$FileName))]
if(length(missing_files) > 0) stop("Missing files:\n", paste(missing_files, collapse="\n"))
message(paste("All", nrow(targets), "files successfully mapped. Proceeding to normalization..."))

# Read and normalize Agilent data
raw_data <- read.maimages(files = targets$FileName, path = txt_path, source = "agilent", green.only = TRUE)
colnames(raw_data) <- paste(targets$Condition, ave(as.character(targets$Condition), targets$Condition, FUN = seq_along), sep="_")

bg_corrected <- backgroundCorrect(raw_data, method="normexp", offset=50)
norm_data <- normalizeBetweenArrays(bg_corrected, method="quantile")

# Filter out control probes
valid_probes <- norm_data$genes$ControlType == 0 
norm_data_filt <- norm_data[valid_probes, ]

# Average replicates
gene_id_col <- ifelse("SystematicName" %in% colnames(norm_data_filt$genes), "SystematicName", "ProbeName")
norm_data_avg <- avereps(norm_data_filt, ID = norm_data_filt$genes[[gene_id_col]])

# Extract the expression matrix for downstream filtering
expr_matrix_final <- norm_data_avg$E

# ==============================================================================
# 3. FILTERING & DIFFERENTIAL EXPRESSION (limma)
# ==============================================================================
# 1. Forgiving human miRNA filter
human_indices <- grep("hsa", rownames(expr_matrix_final), ignore.case = TRUE)
expr_matrix_final <- expr_matrix_final[human_indices, ]

if(nrow(expr_matrix_final) == 0) stop("Error: No human miRNAs (hsa) found. Check annotation mapping.")

# 2. Identify and keep ONLY definitively mapped samples
valid_idx <- which(!is.na(targets$Condition) & targets$Condition %in% c("Normal", "Peritumoral", "Tumor"))

targets <- targets[valid_idx, , drop = FALSE]
expr_matrix_final <- expr_matrix_final[, valid_idx, drop = FALSE]

if(ncol(expr_matrix_final) == 0) stop("Error: No samples mapped to valid conditions.")

# 3. Build the exact-matching design matrix
condition_factor <- factor(targets$Condition, levels=c("Normal", "Peritumoral", "Tumor"))
design <- model.matrix(~0 + condition_factor)
colnames(design) <- c("Normal", "Peritumoral", "Tumor")

# 4. Fit the limma model (3-way contrasts for GSE61710)
fit <- lmFit(expr_matrix_final, design)
cont.matrix <- makeContrasts(TumorVsNormal = Tumor - Normal, 
                             TumorVsPeri = Tumor - Peritumoral, 
                             PeriVsNormal = Peritumoral - Normal, 
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2, trend=TRUE)

# 5. Extract results
res_TvN <- topTable(fit3, coef="TumorVsNormal", adjust="fdr", number=Inf, sort.by="B")
res_TvP <- topTable(fit3, coef="TumorVsPeri", adjust="fdr", number=Inf, sort.by="B")
res_PvN <- topTable(fit3, coef="PeriVsNormal", adjust="fdr", number=Inf, sort.by="B")

write.csv(res_TvN, "GSE61710/GSE61710_Tumor_vs_Normal.csv", row.names = TRUE)
write.csv(res_TvP, "GSE61710/GSE61710_Tumor_vs_Peritumoral.csv", row.names = TRUE)
write.csv(res_PvN, "GSE61710/GSE61710_Peritumoral_vs_Normal.csv", row.names = TRUE)

sig_degs <- res_TvN[res_TvN$adj.P.Val < 0.05 & abs(res_TvN$logFC) > 1, ]
write.csv(sig_degs, "GSE61710/GSE61710_Significant_DEGs_TvN.csv", row.names = TRUE)

# Define unified top 10 miRNAs based on statistical significance (Tumor vs Normal)
top_10_miRNAs <- rownames(res_TvN)[1:10]

# ==============================================================================
# 4. VISUALIZATIONS 
# ==============================================================================
draw_boxplots <- function() {
  par(mfrow=c(1,2), mar=c(10, 4, 4, 2) + 0.1)
  
  # Map dynamic colors based on condition
  color_map <- c("Normal"="blue", "Peritumoral"="orange", "Tumor"="red")
  sample_colors <- color_map[as.character(targets$Condition)]
  
  boxplot(log2(raw_data$E[, valid_idx]), main="Raw Data", las=2, outline=FALSE, col=sample_colors)
  boxplot(expr_matrix_final, main="Normalized Data", las=2, outline=FALSE, col=sample_colors)
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

svg("GSE61710/Plot_01_QC_Boxplots.svg", width=12, height=6)
draw_boxplots()
dev.off()


pca <- prcomp(t(expr_matrix_final), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = targets$Condition)
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition)) +
  geom_point(size=5) + theme_bw() + ggtitle("PCA: GBM vs Normal vs Peritumoral") +
  scale_color_manual(values=c("Normal"="blue", "Peritumoral"="orange", "Tumor"="red"))

ggsave("GSE61710/Plot_02_PCA.svg", plot=pca_plot, width=7, height=6, units="in", dpi=300)

custom_umap <- umap::umap.defaults
custom_umap$n_neighbors <- 5 
umap_out <- umap::umap(t(expr_matrix_final), config=custom_umap)
umap_df <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Condition = targets$Condition)
umap_plot <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=Condition)) +
  geom_point(size=5) + theme_bw() + ggtitle("UMAP: GBM vs Normal vs Peritumoral") +
  scale_color_manual(values=c("Normal"="blue", "Peritumoral"="orange", "Tumor"="red"))
ggsave("GSE61710/Plot_03_UMAP.svg", plot=umap_plot, width=7, height=6, units="in", dpi=300)

svg("GSE61710/Plot_04_Mean_Variance.svg", width=8, height=6)
plotSA(fit3, main="Mean-Variance Trend", ylab="sqrt(sigma)", xlab="Average log-expression")
dev.off()

svg("GSE61710/Plot_05_QQ_Plot.svg", width=8, height=6)
qqt(fit3$t, df = fit3$df.total + fit3$df.prior, main="Q-Q Plot", ylab="Sample Quantiles", xlab="Theoretical Quantiles")
abline(0,1, col="red")
dev.off()

svg("GSE61710/Plot_06_MD_Plot.svg", width=8, height=6)
plotMD(fit3, column=1, status=decideTests(fit3), main="MD Plot: GBM vs Normal")
dev.off()

volc <- EnhancedVolcano(res_TvN, 
                        lab = rownames(res_TvN), 
                        x = 'logFC', 
                        y = 'adj.P.Val', 
                        selectLab = top_10_miRNAs,
                        title = 'GSE61710: GBM vs Normal', 
                        pCutoff = 0.05, 
                        FCcutoff = 1.0, 
                        pointSize = 3.0, 
                        labSize = 3.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5)
ggsave("GSE61710/Plot_07_Volcano.svg", plot=volc, width=12, height=10, units="in", dpi=300)

svg("GSE61710/Plot_08_Adj_PVal_Hist.svg", width=8, height=6)
hist(res_TvN$adj.P.Val, col="lightblue", main="Histogram of Adj P-values", xlab="Adjusted P-value")
dev.off()

results_decision <- decideTests(fit3, adjust.method="fdr", p.value=0.05)
colnames(results_decision) <- gsub("Tumor", "GBM", colnames(results_decision))
svg("GSE61710/Plot_09_VennDiagram.svg", width=8, height=8)
vennDiagram(results_decision, include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("red", "orange", "blue"))
dev.off()

top_50 <- rownames(res_TvN)[1:min(50, nrow(res_TvN))]
heatmap_mat <- expr_matrix_final[top_50, ]
svg("GSE61710/Plot_10_Heatmap.svg", width=8, height=10)
pheatmap(heatmap_mat, scale = "row", main = "Top 50 Differentially Expressed miRNAs", annotation_col = data.frame(Condition = targets$Condition, row.names=colnames(heatmap_mat)), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

svg("GSE61710/Plot_11_Density_Plot.svg", width=12, height=6)
par(mfrow=c(1,2))
# Map dynamic colors for densities
color_map <- c("Normal"="blue", "Peritumoral"="orange", "Tumor"="red")
density_colors <- color_map[unique(as.character(targets$Condition))]

plotDensities(log2(raw_data$E[, valid_idx]), group = targets$Condition, main = "Raw Density", legend = "topright", col = density_colors)
plotDensities(expr_matrix_final, group = targets$Condition, main = "Normalized Density", legend = "topright", col = density_colors)
par(mfrow=c(1,1))
dev.off()

normal_idx <- which(targets$Condition == "Normal")
tumor_idx <- which(targets$Condition == "Tumor")
mean_normal <- if(length(normal_idx) > 1) rowMeans(expr_matrix_final[, normal_idx]) else expr_matrix_final[, normal_idx]
mean_tumor <- if(length(tumor_idx) > 1) rowMeans(expr_matrix_final[, tumor_idx]) else expr_matrix_final[, tumor_idx]
sig_genes <- rownames(res_TvN)[res_TvN$adj.P.Val < 0.05 & abs(res_TvN$logFC) > 1]
is_sig <- rownames(expr_matrix_final) %in% sig_genes

svg("GSE61710/Plot_12_Scatter_Means.svg", width=8, height=8)
plot(mean_normal, mean_tumor, pch = 16, cex = 0.5, col = "grey", 
     xlab = "Average Normal Expression (Log2)", 
     ylab = "Average GBM Expression (Log2)", 
     main = "Scatter Plot: Mean Normal vs Mean GBM")
abline(0, 1, col = "blue", lwd = 2)
points(mean_normal[is_sig], mean_tumor[is_sig], col = "red", pch = 16, cex = 0.6)
top_10_idx <- match(top_10_miRNAs, rownames(expr_matrix_final))
text(mean_normal[top_10_idx], mean_tumor[top_10_idx], labels = top_10_miRNAs, pos = 4, cex = 0.8, font = 2, col = "black")
dev.off()

draw_scatter <- function() {
  scatter_df <- data.frame(
    miRNA = as.character(rownames(expr_matrix_final)),
    Mean_Normal = as.numeric(mean_normal),
    Mean_GBM = as.numeric(mean_tumor),
    Significance = ifelse(is_sig, "Significant", "Not Significant"),
    stringsAsFactors = FALSE
  )
  
  # Use NA so ggrepel cleanly ignores non-target miRNAs
  scatter_df$Label <- ifelse(scatter_df$miRNA %in% top_10_miRNAs, scatter_df$miRNA, NA)
  
  p <- ggplot(scatter_df, aes(x = Mean_Normal, y = Mean_GBM)) +
    geom_point(aes(color = Significance), size = 1.5, alpha = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 1) +
    geom_text_repel(aes(label = Label),
                    na.rm = TRUE, # Drops all NA rows without errors
                    size = 3.5,
                    fontface = "bold",
                    box.padding = 0.8,
                    point.padding = 0.3,
                    segment.color = "black",
                    min.segment.length = 0, 
                    max.overlaps = Inf) +
    theme_bw() +
    labs(x = "Average Normal Expression (Log2)",
         y = "Average GBM Expression (Log2)",
         title = "Scatter Plot: Mean Normal vs Mean GBM") +
    theme(legend.position = "none")
  
  return(p)
}

scatter_plot <- draw_scatter()
ggsave("GSE61710/Plot_12_Scatter_Means_Connector.svg", plot = scatter_plot, width = 8, height = 8,units = "in")

# ==============================================================================
# 5. MACHINE LEARNING DATA PREPARATION (Normal & Tumor ONLY)
# ==============================================================================
ml_idx <- which(targets$Condition %in% c("Normal", "Tumor"))
expr_ml <- expr_matrix_final[, ml_idx, drop = FALSE]
targets_ml <- targets[ml_idx, , drop = FALSE]
ml_matrix <- t(expr_ml)
z_score_matrix <- scale(ml_matrix)

ml_df <- as.data.frame(z_score_matrix, check.names = FALSE)
ml_df$Target <- ifelse(targets_ml$Condition == "Tumor", 1, 0)
ml_df$Sample_ID <- rownames(ml_df)
ml_df <- ml_df[, c("Sample_ID", colnames(z_score_matrix), "Target")]
print(table(ml_df$Target))
write.csv(ml_df, "GSE61710/GSE61710_ML_Ready.csv", row.names = FALSE)