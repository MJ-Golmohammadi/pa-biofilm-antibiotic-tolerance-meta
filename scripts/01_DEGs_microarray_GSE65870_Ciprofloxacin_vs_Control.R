# create new folder for results and figures
dir.create("results/GSE65870", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/GSE65870", showWarnings = FALSE, recursive = TRUE)

# Clear the workspace
rm(list = ls())

library(GEOquery)
library(affy)
library(limma)

# Unpack the raw data
untar(
  "data/GSE65870_RAW.tar",
  exdir = "data/GSE65870"
)


# Load raw microarray data
raw_data <- ReadAffy(
  celfile.path = "data/GSE65870"
)
colnames(raw_data)
# Normalization using RMA
eset <- rma(raw_data)
expr <- exprs(eset)
colnames(expr)
summary(expr)
#################################################
# Filter low-expression probes
#################################################

keep <- rowSums(expr > 3) >= 2

expr <- expr[
  keep,
]

cat(
  "Remaining probes:",
  nrow(expr),
  "\n"
)
# summary of expr
summary(expr)


# Define the group labels for the samples
group <- factor(c(
  rep("Ciprofloxacin_Treatment",3),
  rep("Tobramycin_Treatment",3),
  rep("Control",3)
))

# Check the dimensions of the expression matrix
dim(expr)


################################################## 
# Boxplot of normalized expression data
#################################################

# Set column names to sample names without .CEL.gz extension
colnames(expr) <- sub("\\.CEL\\.gz$", "", colnames(expr))

# Create a boxplot of the normalized expression data
#################################################
# Colors by group
#################################################

group_colors <- c(
  rep("#56B4E9",3),   # Ciprofloxacin
  rep("#E69F00",3),   # Tobramycin
  rep("#CC79A7",3)    # Control
)

#################################################
# Boxplot
#################################################

png(
  "figures/GSE65870/RMA_normalized_expression_boxplot_GSE65870.png",
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  expr,
  las = 2,
  col = group_colors,
  border = "grey30",
  outline = FALSE,
  main = "RMA-normalized expression distribution (GSE65870)",
  ylab = "Log2 expression",
  cex.axis = 0.9,
  cex.main = 1.2
)

legend(
  #add position for legend outside the plot
  "topright",
  # add legend size
  cex = 0.8, 
  legend = c(
    "Ciprofloxacin_Treatment",
    "Tobramycin_Treatment",
    "Control"
  ),
  fill = c(
    "#56B4E9",
    "#E69F00",
    "#CC79A7"
  )
)

dev.off()

# Save the normalized expression data to a CSV file
write.csv(
  expr,
  "results/GSE65870/GSE65870_RMA_expression.csv"
)

colnames(expr)
length(group)
ncol(expr)
#################################################
# PCA 
#################################################

library(ggplot2)

pca <- prcomp(t(expr), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  condition = group
)

percentVar <- pca$sdev^2 / sum(pca$sdev^2)

p <- ggplot(pca_df,
            aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  labs(
    title = "RMA normalized expression PCA - GSE65870",
    x = paste0("PC1: ", round(percentVar[1] * 100, 1), "% variance"),
    y = paste0("PC2: ", round(percentVar[2] * 100, 1), "% variance")
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p
#save PCA plot to figures folder
ggsave(
  "figures/GSE65870/PCA_GSE65870.png",
  plot = p,
  width = 8,
  height = 6
)

#################################################
## Sample distance heatmap
#################################################

library(pheatmap)

# Calculate sample-to-sample distances
sampleDists <- dist(t(expr))

# Convert distance object to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Sample annotation
annotation_col <- data.frame(
  Condition = factor(group)
)

rownames(annotation_col) <- colnames(expr)

# Draw and save heatmap
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  main = "Sample-to-Sample Distance - GSE65870",
  filename = "figures/GSE65870/Sample_distance_heatmap_GSE65870.png",
  width = 6,
  height = 5,
  fontsize = 8,       
  fontsize_row = 7,    
  fontsize_col = 7 
)


# Create the design matrix for differential expression analysis
design <- model.matrix(~0 + group)
# Name the columns of the design matrix according to the group levels
colnames(design) <- levels(group)

design

# Define the contrast matrix for the comparison of interest
contrast.matrix <- makeContrasts(
  Ciprofloxacin_vs_Control =
    Ciprofloxacin_Treatment - Control,

  levels = design
)


library(limma)
# Fit the linear model and perform the differential expression analysis
fit <- lmFit(expr, design)
# Apply the contrast matrix to the fitted model
fit2 <- contrasts.fit(fit, contrast.matrix)
# Apply empirical Bayes moderation to the fitted model
fit2 <- eBayes(fit2)
# Extract the top differentially expressed genes based on the specified contrast
deg <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH"
)

head(deg)
dim(deg)

deg_sorted <- deg[order(deg$adj.P.Val), ]

write.csv(
  deg_sorted,
  "results/GSE65870/GSE65870_DEGs_Ciprofloxacin_vs_Control.csv"
)


head(rownames(deg))


#################################################
# Annotation
#################################################

library(readr)

annot <- read_csv(
  "data/Pseudomonas_aeruginosa_PAO1_107.csv",
  show_col_types = FALSE
)

dim(annot)
colnames(annot)

############################################################
# Keep only first 6 characters of probe IDs
############################################################

deg$Locus_Tag <- substr(
  rownames(deg),
  1,
  6
)


head(
  data.frame(
    Original = rownames(deg),
    ID = deg$Locus_Tag
  )
)


deg_annot <- merge(
  deg,
  annot,
  by = "Locus_Tag",
  all.x = TRUE
)


sum(is.na(deg_annot$Gene_Name))
sum(is.na(deg_annot$Locus_Tag))


deg_final <- deg_annot[, c(
  "Locus_Tag",
  "Gene_Name",
  "logFC",
  "AveExpr",
  "P.Value",
  "adj.P.Val",
  "B"
)]

deg_final <- deg_final[order(
    deg_annot$adj.P.Val,
    -abs(deg_annot$logFC)
), ]

write.csv(
  deg_final,
  "results/GSE65870/GSE65870_DEGs_Annotated_Ciprofloxacin_vs_Control.csv",
  row.names = FALSE
)

deg_final[1:20,
          c("Locus_Tag",
            "Gene_Name",
            "logFC",
            "adj.P.Val")]


#################################################
# Create final Meta-analysis input file
#################################################

# Gene column:
# if Gene_Name is empty use Locus_Tag

deg_final$Gene <- ifelse(
  is.na(deg_final$Gene_Name) |
    deg_final$Gene_Name == "",
  deg_final$Locus_Tag,
  deg_final$Gene_Name
)

#################################################
# Keep important columns
#################################################

meta_deg <- deg_final[, c(
  "Gene",
  "Locus_Tag",
  "logFC",
  "adj.P.Val"
)]


#################################################
# Sort by adjusted p-value and log fold change
#################################################

meta_deg <- meta_deg[
  order(meta_deg$adj.P.Val,
        -abs(meta_deg$logFC)),
]

head(meta_deg)
#################################################
# Remove duplicated genes
#################################################

meta_deg <- meta_deg[!duplicated(meta_deg$Gene), ]

#################################################
# Sort by adjusted p-value
#################################################

meta_deg <- meta_deg[
  order(meta_deg$adj.P.Val),
] 

#################################################
# Save file
#################################################

write.csv(
  meta_deg,
  "results/GSE65870/GSE65870_meta_input_Ciprofloxacin_vs_Control.csv",
  row.names = FALSE
)

#################################################
# Significant File
#################################################

# Filter the DEGs based on adjusted p-value and log fold change thresholds
deg_sig <- subset(
  meta_deg,
  adj.P.Val < 0.05 & abs(logFC) > 1
)


write.csv(deg_sig, "results/GSE65870/GSE65870_significant_Ciprofloxacin_vs_Control.csv")
#################################################
# Statistics
#################################################

cat("Total genes:", nrow(meta_deg), "\n")

cat(
  "Significant genes:",
  sum(
    meta_deg$adj.P.Val < 0.05 &
      abs(meta_deg$logFC) > 1
  ),
  "\n"
)


#################################################
# Upregulated
#################################################

up_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC > 1
)

#################################################
# Downregulated
#################################################

down_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC < -1
)

cat("Upregulated:", nrow(up_genes), "\n")
cat("Downregulated:", nrow(down_genes), "\n")

#################################################
# Save significant genes
#################################################

write.csv(
  up_genes,
  "results/GSE65870/GSE65870_upregulated_Ciprofloxacin_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  down_genes,
  "results/GSE65870/GSE65870_downregulated_Ciprofloxacin_vs_Control.csv",
  row.names = FALSE
)


#################################################
# Volcano Plot
#################################################

library(EnhancedVolcano)
library(ggplot2)

topGenes <- head(
  meta_deg$Gene,
  10
)


(
  EnhancedVolcano(
      meta_deg,
  lab = meta_deg$Gene,
  selectLab = topGenes,
  x = "logFC",
  y = "adj.P.Val",

  title =
    "Ciprofloxacin-treated biofilm vs Control biofilm (GSE65870)",

  subtitle =
    "",

  xlab = expression(Log[2]~Fold~Change),
  ylab = expression(-Log[10]~Adjusted~P),

    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 6,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1,
    colConnectors = 'black',
    max.overlaps = 20,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 1.2,
    borderColour = 'black'
  )
)+
    coord_cartesian(ylim = c(0, 8.5)) +
      theme(
      plot.title = element_text(hjust = 0.5),     
      plot.subtitle = element_text(hjust = 0.5) 
      )


ggsave(
  "figures/GSE65870/Volcano_plot_GSE65870_Ciprofloxacin_vs_Control.png",
  width = 9.5,
  height = 12,
  dpi = 300
)

#################################################
# MA Plot - Publication Quality
#################################################

library(ggplot2)
library(ggrepel)

#################################################
# Prepare data
#################################################

ma_df <- deg_final

#################################################
# Sort by adjusted p-value and log fold change
#################################################

deg_final <- deg_final[
  order(deg_final$adj.P.Val,
        -abs(deg_final$logFC)),
]

#################################################
# Remove duplicated genes
#################################################

deg_final <- deg_final[!duplicated(deg_final$Gene), ]
head(deg_final)

#################################################
# Classification
#################################################

ma_df$Status <- "Not Significant"

ma_df$Status[
  ma_df$adj.P.Val < 0.05 &
    ma_df$logFC > 1
] <- "Upregulated"

ma_df$Status[
  ma_df$adj.P.Val < 0.05 &
    ma_df$logFC < -1
] <- "Downregulated"

#################################################
# Top genes for labeling
#################################################

top_genes <- ma_df[
  order(ma_df$adj.P.Val),
][1:10, ]

#################################################
# Plot
#################################################

p_ma <- ggplot(
  ma_df,
  aes(
    x = AveExpr,
    y = logFC
  )
) +

  geom_point(
    aes(color = Status),
    alpha = 0.75,
    size = 2
  ) +

  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +

  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2
  ) +

  scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Not Significant" = "grey80"
    )
  ) +

  labs(
    title =
      "Ciprofloxacin-treated biofilm vs Control biofilm (GSE65870)",

    subtitle =
      "",

    x = "Average Expression (A)",

    y = expression(Log[2]~Fold~Change),

    color = NULL
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    legend.position = "top",

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      linewidth = 0.8
    )
  )

#################################################
# Show plot
#################################################

p_ma

#################################################
# Save figure
#################################################

ggsave(
  "figures/GSE65870/MA_plot_GSE65870_Ciprofloxacin_vs_Control.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 300
)
