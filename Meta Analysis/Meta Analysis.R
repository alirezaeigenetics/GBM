# ==============================================================================
# 1. SETUP AND PACKAGE LOADING
# ==============================================================================
# Install required packages if missing
packages_needed <- c("ggVennDiagram", "UpSetR", "ggplot2", "dplyr", "ggrepel")
for (p in packages_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Use dynamic path resolution anchoring to the R Project root directory
script_dir <- getwd()
cat("\nRunning Meta-Analysis from root directory:", script_dir, "\n")

# ==============================================================================
# 2. LOAD INDIVIDUAL DATASET CSVs FROM THEIR FOLDERS
# ==============================================================================
DEG_M1 <- read.csv("GSE61710/GSE61710_Tumor_vs_Normal.csv", row.names = 1)
DEG_M2 <- read.csv("GSE65626/GSE65626_GBM_vs_Normal.csv", row.names = 1)
DEG_M3 <- read.csv("GSE90603/GSE90603_GBM_vs_Normal.csv", row.names = 1)
DEG_M4 <- read.csv("GSE91014/GSE91014_Tumor_vs_Normal.csv", row.names = 1)
DEG_M5 <- read.csv("GSE93850/GSE93850_Tumor_vs_Normal.csv", row.names = 1)
DEG_M6 <- read.csv("GSE103228/GSE103228_Tumor_vs_Normal.csv", row.names = 1)
DEG_TCGA <- read.csv("TCGA_RNAseq/TCGA_DE_miRNAs.csv", row.names = 1)

# ==============================================================================
# 3. DEFINE SMART PROCESSING FUNCTIONS (NO FILTERING)
# ==============================================================================
process_microarray <- function(df) {
  mirnas <- rownames(df)
  # Normalize names to match TCGA format: lowercase, remove -3p/-5p/*
  mirnas_norm <- gsub("-[35]p$|\\*$", "", tolower(mirnas))
  return(unique(mirnas_norm))
}

process_tcga <- function(df) {
  mirnas <- rownames(df)
  # Normalize names to lowercase and replace dots with hyphens
  mirnas_norm <- gsub("\\.", "-", tolower(mirnas))
  return(unique(mirnas_norm))
}

# Apply functions to extract & normalize miRNAs for UpSet/Venn
all_M1 <- process_microarray(DEG_M1)
all_M2 <- process_microarray(DEG_M2)
all_M3 <- process_microarray(DEG_M3)
all_M4 <- process_microarray(DEG_M4)
all_M5 <- process_microarray(DEG_M5)
all_M6 <- process_microarray(DEG_M6)
all_TCGA <- process_tcga(DEG_TCGA)

# Bundle them into a list for the plots
All_list <- list(
  GSE61710 = all_M1,
  GSE65626 = all_M2,
  GSE90603 = all_M3,
  GSE91014 = all_M4,
  GSE93850 = all_M5,
  GSE103228 = all_M6,
  TCGA_GBM = all_TCGA
)

# ==============================================================================
# 4. OVERLAP VISUALIZATIONS (Venn & UpSet - PNG & SVG)
# ==============================================================================
# Ensure Meta Analysis folder exists for outputs
meta_dir <- "Meta Analysis"
if(!dir.exists(meta_dir)) dir.create(meta_dir)

# The Venn Diagram
venn_plot <- ggVennDiagram(All_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#E7B800") +
  theme(legend.position = "none") +
  ggtitle("Overlap of ALL Detected miRNAs across 7 GBM Datasets")

ggsave(file.path(meta_dir, "MetaAnalysis_All_VennDiagram.png"), plot = venn_plot, width = 12, height = 10, units = "in", dpi = 300, bg = "white")
ggsave(file.path(meta_dir, "MetaAnalysis_All_VennDiagram.svg"), plot = venn_plot, width = 12, height = 10, units = "in", bg = "white")

# The UpSet Plot
upset_plot <- upset(fromList(All_list), 
                    order.by = "freq", 
                    nsets = 7, 
                    nintersects = 20, 
                    main.bar.color = "#2E9FDF", 
                    sets.bar.color = "#E7B800",
                    text.scale = c(1.5, 1.2, 1.5, 1.2, 1.2, 1.5))

png(file.path(meta_dir, "MetaAnalysis_All_UpSetPlot.png"), width = 10, height = 6, units = "in", res = 300)
print(upset_plot)
dev.off()

svg(file.path(meta_dir, "MetaAnalysis_All_UpSetPlot.svg"), width = 10, height = 6)
print(upset_plot)
dev.off()

# Extract universally detected miRNAs (present in all 7 datasets)
common_all_mirnas <- Reduce(intersect, All_list)
cat("\nTotal miRNAs detected in all 7 datasets:", length(common_all_mirnas), "\n")
write.csv(common_all_mirnas, file.path(meta_dir, "Conserved_All_miRNAs.csv"), row.names = FALSE)

# ==============================================================================
# 5. STANDARDIZATION ENGINE FOR FISHER'S META-ANALYSIS
# ==============================================================================
standardize_deg <- function(df, is_tcga = FALSE) {
  if(is_tcga) {
    df$miRNA <- gsub("\\.", "-", tolower(rownames(df)))
    df$logFC <- df$log2FoldChange
    df$P_Value <- df$padj
  } else {
    df$miRNA <- gsub("-[35]p$|\\*$", "", tolower(rownames(df)))
    df$logFC <- df$logFC
    df$P_Value <- df$adj.P.Val
  }
  
  df <- df[!is.na(df$P_Value), ]
  return(df[, c("miRNA", "logFC", "P_Value")])
}

# Standardize all datasets
deg_list <- list(
  GSE61710 = standardize_deg(DEG_M1, FALSE), 
  GSE65626 = standardize_deg(DEG_M2, FALSE), 
  GSE90603 = standardize_deg(DEG_M3, FALSE), 
  GSE91014 = standardize_deg(DEG_M4, FALSE), 
  GSE93850 = standardize_deg(DEG_M5, FALSE), 
  GSE103228 = standardize_deg(DEG_M6, FALSE), 
  TCGA = standardize_deg(DEG_TCGA, TRUE)
)

# ==============================================================================
# 6. RUN FISHER'S META-ANALYSIS 
# ==============================================================================
cat("\nRunning Custom Fisher's Meta-Analysis...\n")
all_deg <- bind_rows(deg_list, .id = "Study")

meta_results <- all_deg %>%
  group_by(miRNA) %>%
  summarise(
    Studies_Found = n(),
    Avg_logFC = mean(logFC, na.rm = TRUE),
    Chi_Square = -2 * sum(log(P_Value), na.rm = TRUE),
    Meta_P_Value = pchisq(Chi_Square, df = 2 * n(), lower.tail = FALSE)
  )

final_meta_mirnas <- meta_results %>%
  filter(Meta_P_Value < 0.05 & Studies_Found >= 4) %>%
  arrange(Meta_P_Value)

write.csv(final_meta_mirnas, file.path(meta_dir, "MetaAnalysis_Fisher_Results.csv"), row.names = FALSE)

# ==============================================================================
# 7. CUSTOM META-VOLCANO PLOT (PNG & SVG)
# ==============================================================================
meta_results$Significance <- ifelse(meta_results$Meta_P_Value < 0.05 & meta_results$Avg_logFC > 0.5, "Upregulated",
                                    ifelse(meta_results$Meta_P_Value < 0.05 & meta_results$Avg_logFC < -0.5, "Downregulated", "Not Significant"))

# Target miRNAs for labeling based on final ML model
target_mirnas <- c("hsa-mir-137", "hsa-mir-765", "hsa-mir-485", "hsa-mir-769")
watch_list_data <- meta_results %>% filter(miRNA %in% target_mirnas)

volcano_meta_plot <- ggplot(meta_results, aes(x = Avg_logFC, y = -log10(Meta_P_Value), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Downregulated" = "#2E9FDF", "Not Significant" = "grey", "Upregulated" = "#E7B800")) +
  geom_label_repel(data = watch_list_data, 
                   aes(label = miRNA), 
                   color = "black",       
                   box.padding = 0.5,     
                   point.padding = 0.5,   
                   show.legend = FALSE) + 
  theme_minimal() +
  labs(title = "Fisher's Meta-Analysis Volcano Plot",
       subtitle = "Combined P-values across 7 GBM Datasets",
       x = "Average log2 Fold Change",
       y = "-log10(Meta P-Value)") +
  theme(text = element_text(size = 14),
        legend.position = "bottom") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

ggsave(file.path(meta_dir, "MetaAnalysis_Custom_Volcano_Labeled.png"), plot = volcano_meta_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
ggsave(file.path(meta_dir, "MetaAnalysis_Custom_Volcano_Labeled.svg"), plot = volcano_meta_plot, width = 10, height = 8, units = "in", bg = "white")

cat("\nMeta-Analysis Complete! PNG and SVG outputs saved in 'Meta Analysis' folder.\n")