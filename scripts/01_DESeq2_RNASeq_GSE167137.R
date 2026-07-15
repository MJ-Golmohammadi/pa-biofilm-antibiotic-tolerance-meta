# Clear the workspace
rm(list = ls())

library(DESeq2)
library(readxl)
library(dplyr)
library(readr)

#################################################
# Read raw counts sheet
#################################################

counts <- read_xlsx(
  "data/GSE167137/GSE167137_P_aeruginosa_count_data.xlsx",
  sheet = "Raw counts"
)
colnames(counts)

############################################################
# Read PAO1 annotation
############################################################

annot <- read_csv(
  "data/Pseudomonas_aeruginosa_PAO1_107.csv",
  show_col_types = FALSE
)

############################################################
# Keep only required annotation columns
############################################################

annot_map <- annot %>%
  select(
    Locus_Tag,
    Gene_Name
  ) %>%
  distinct()

############################################################
# Rename gene ID column
############################################################

colnames(counts)[1] <- "Gene ID"
colnames(counts)

############################################################
# Replace Gene IDs with Gene Names when available
############################################################

counts_new <- counts %>%
  left_join(
    annot_map,
    by = c("Gene ID" = "Locus_Tag")
  ) %>%
  mutate(
    `Gene ID` = ifelse(
      !is.na(Gene_Name) & Gene_Name != "",
      Gene_Name,
      `Gene ID`
    )
  ) %>%
  select(
    -Gene_Name
  )


############################################################
# Summary
############################################################

cat(
  "Mapped genes:",
  sum(!is.na(
    left_join(
      counts,
      annot_map,
      by = c("Gene ID" = "Locus_Tag")
    )$Gene_Name
  )),
  "\n"
)

############################################################
# Save updated file
############################################################

write.csv(
  counts_new,
  "data/GSE167137/GSE167137_count_data_GeneSymbol.csv",
  row.names = FALSE
)


################################################
# convert counts to counts_new
################################################
counts <- counts_new
sum(duplicated(counts$`Gene ID`))

############################################################
# Collapse duplicated gene symbols
############################################################

library(dplyr)

counts <- counts %>%
  group_by(`Gene ID`) %>%
  summarise(
    across(
      where(is.numeric),
      sum,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

############################################################
# Set gene IDs as row names
############################################################

rownames(counts) <- counts$`Gene ID`

############################################################
# Remove Gene ID column for DESeq2
############################################################

counts_matrix <- counts %>%
  select(-`Gene ID`)

############################################################
# Quality control
############################################################

cat(
  "Duplicated genes:",
  sum(duplicated(rownames(counts_matrix))),
  "\n"
)

#################################################
# Keep PA samples only
#################################################

counts_pa <- counts[, c(
  "Gene ID",
  "PA-0M-4",
  "PA-0M-6",
  "PA-0M-7",
  "PA-0M-8",
  "PA-5M-4",
  "PA-5M-6",
  "PA-5M-7",
  "PA-5M-8"
)]

#################################################
# Convert to matrix
#################################################

counts_pa <- as.data.frame(counts_pa)

rownames(counts_pa) <- counts_pa$`Gene ID`
counts_pa$`Gene ID` <- NULL

count_matrix <- as.matrix(counts_pa)

storage.mode(count_matrix) <- "integer"

keep <- rowSums(count_matrix >= 10) >= 2
count_matrix <- count_matrix[keep, ]
summary(count_matrix)
#################################################
# Metadata
#################################################

condition <- factor(c(
  rep("Control",4),
  rep("Meropenem",4)
))

coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition
)

#################################################
# DESeq2
#################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

#################################################
# VST Boxplot - Publication Quality
#################################################

library(reshape2)
library(ggplot2)

#################################################
# Extract VST matrix
#################################################

vst_mat <- assay(vsd)

#################################################
# Convert to data.frame
#################################################

vst_df <- as.data.frame(vst_mat)

vst_df$Gene <- rownames(vst_df)

#################################################
# Long format
#################################################

vst_long <- melt(
  vst_df,
  id.vars = "Gene",
  variable.name = "Sample",
  value.name = "Expression"
)

#################################################
# Sample annotation
#################################################

vst_long$Condition <- ifelse(
  grepl("^PA-0M", vst_long$Sample),
  "Control",
  "Meropenem"
)

vst_long$Condition <- factor(
  vst_long$Condition,
  levels = c(
    "Control",
    "Meropenem"
  )
)

#################################################
# Boxplot
#################################################

p_box <- ggplot(
  vst_long,
  aes(
    x = Sample,
    y = Expression,
    fill = Condition
  )
) +

  geom_boxplot(
    linewidth = 0.4,
    outlier.size = 0.5
  ) +

  scale_fill_manual(
    values = c(
      "Control" = "#56B4E9",
      "Meropenem" = "#E69F00"
    )
  ) +

  labs(
    title = "VST-normalized expression distribution (GSE167137)",
    subtitle = "DESeq2 variance stabilizing transformation",
    x = NULL,
    y = "VST expression"
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),

    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),

    axis.text.y = element_text(
      size = 10
    ),

    legend.position = "top",

    legend.title = element_blank(),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      linewidth = 0.8
    )
  )

#################################################
# Show plot
#################################################

p_box

#################################################
# Save figure
#################################################

ggsave(
  "figures/GSE167137/VST_boxplot_GSE167137.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 300
)


#################################################
# PCA
#################################################

# create new folder for results
dir.create("results/GSE167137", showWarnings = FALSE)

# create new folder for figures
dir.create("figures/GSE167137", showWarnings = FALSE)

# Save the normalized expression data to a CSV file
write.csv(assay(vsd), "results/GSE167137/VST_expression.csv")

p <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA of VST-normalized counts - GSE167137") +
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
  "figures/GSE167137/PCA_plot_GSE167137.png",
  width = 8,
  height = 6
)

#################################################
## Sample distance heatmap
#################################################

library(pheatmap)

# Sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Sample annotation
annotation_col <- data.frame(
  Condition = colData(vsd)$condition
)

rownames(annotation_col) <- colnames(vsd)

# Draw and save heatmap
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-Sample Distance after VST normalization - GSE167137",
  filename = "figures/GSE167137/Sample_distance_heatmap_GSE167137.png",
  width = 6,
  height = 5,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Differential expression
#################################################

res <- results(
  dds,
  contrast = c(
    "condition",
    "Meropenem",
    "Control"
  )
)

#################################################
# Summary
#################################################

resultsNames(dds)

summary(res)



#################################################
# Convert results to data frame
#################################################

deg <- as.data.frame(res)

deg$Gene <- rownames(deg)


#################################################
# Remove NA adjusted p-values
#################################################

deg <- deg[
  !is.na(deg$padj),
]

#################################################
# Sort by adjusted p-value
#################################################

deg <- deg[
  order(deg$padj),
]

#################################################
# remove duplicates (keep best one)
#################################################

deg <- deg[!duplicated(deg$Gene), ]


#################################################
# Sort by adjusted p-value
#################################################

deg <- deg[
  order(deg$padj),
]

#################################################
# Save full DEG table
#################################################

write.csv(
  deg,
  "results/GSE167137/GSE167137_DEGs.csv",
  row.names = FALSE
)


#################################################
# Create meta-analysis input file
#################################################

meta_deg <- deg[, c(
  "Gene",
  "log2FoldChange",
  "padj"
)]

colnames(meta_deg) <- c(
  "Gene",
  "logFC",
  "adj.P.Val"
)

#################################################
# Check for duplicates
#################################################

sum(duplicated(rownames(counts)))
sum(duplicated(meta_deg$Gene))

#################################################
# Sort by adjusted p-value and log fold change
#################################################

meta_deg <- meta_deg[
  order(meta_deg$adj.P.Val,
        -abs(meta_deg$logFC)),
]

#################################################
# Save meta-analysis file
#################################################

write.csv(
  meta_deg,
  "results/GSE167137/GSE167137_meta_input.csv",
  row.names = FALSE
)

#################################################
# Significant DEGs
#################################################

deg_sig <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

cat(
  "Significant DEGs:",
  nrow(deg_sig),
  "\n"
)
cat(
  "Upregulated:",
  nrow(subset(deg_sig, logFC > 1)),
  "\n"
)
cat(
  "Downregulated:",
  nrow(subset(deg_sig, logFC < -1)),
  "\n"
)

write.csv(
  deg_sig,
  "results/GSE167137/GSE167137_significant.csv",
  row.names = FALSE
)
#################################################
# Summary statistics
#################################################

cat(
  "Total genes:",
  nrow(meta_deg),
  "\n"
)

cat(
  "Significant genes:",
  sum(
    meta_deg$adj.P.Val < 0.05 &
      abs(meta_deg$logFC) > 1,
    na.rm = TRUE
  ),
  "\n"
)

sum(
  meta_deg$adj.P.Val < 0.05,
  na.rm = TRUE
)
resultsNames(dds)

summary(res)

head(
  deg[,c(
    "Gene",
    "log2FoldChange",
    "pvalue",
    "padj"
  )],
  20
)
#################################################
# Upregulated genes
#################################################

up_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC > 1
)

#################################################
# Downregulated genes
#################################################

down_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC < -1
)

cat(
  "Upregulated:",
  nrow(up_genes),
  "\n"
)

cat(
  "Downregulated:",
  nrow(down_genes),
  "\n"
)

#################################################
# Save upregulated genes
#################################################

write.csv(
  up_genes,
  "results/GSE167137/GSE167137_upregulated.csv",
  row.names = FALSE
)

#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE167137/GSE167137_downregulated.csv",
  row.names = FALSE
)

#################################################
# Show top genes
#################################################

head(meta_deg, 20)

#################################################
# DESeq2 summary
#################################################

summary(res)


#################################################
# LFC Shrinkage
#################################################

library(apeglm)

res_shrink <- lfcShrink(
  dds,
  coef = "condition_Meropenem_vs_Control",
  type = "apeglm"
)

#################################################
# Convert to data frame
#################################################

deg_shrink <- as.data.frame(res_shrink)

deg_shrink$Gene <- rownames(deg_shrink)

deg_shrink <- deg_shrink[
  !is.na(deg_shrink$padj),
]

deg_shrink <- deg_shrink[
  order(deg_shrink$padj),
]

#################################################
# Top genes for labels
#################################################

topGenes <- head(
  deg_shrink$Gene,
  10
)

#################################################
# Volcano Plot
#################################################

library(EnhancedVolcano)
library(ggplot2)


(
  EnhancedVolcano(
    deg_shrink,
    lab = deg_shrink$Gene,
    selectLab = topGenes, 
    x = 'log2FoldChange',
    y = 'padj',
    xlab = 'log2 Fold Change',
    ylab = '-log10 Adjusted P-value',
    title = 'Meropenem-treated biofilm vs Control biofilm (GSE167137)',
    subtitle = '',
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
) +
    theme(
    plot.title = element_text(hjust = 0.5),     
    plot.subtitle = element_text(hjust = 0.5) 
    )

#################################################
# Save PNG
#################################################

ggsave(
  "figures/GSE167137/Volcano_plot_shrink_GSE167137.png",
  width = 9.5,
  height = 9.5,
  dpi = 300
)

#################################################
# MA Plot - Publication Quality
#################################################

library(ggplot2)
library(ggrepel)

# Convert to data frame
ma_df <- as.data.frame(res_shrink)

ma_df$Gene <- rownames(ma_df)

# Remove NAs
ma_df <- ma_df[
  !is.na(ma_df$padj) &
  !is.na(ma_df$baseMean),
]

#################################################
# Classification
#################################################

ma_df$Status <- "Not Significant"

ma_df$Status[
  ma_df$padj < 0.05 &
    ma_df$log2FoldChange > 1
] <- "Upregulated"

ma_df$Status[
  ma_df$padj < 0.05 &
    ma_df$log2FoldChange < -1
] <- "Downregulated"

#################################################
# Top genes for labeling
#################################################

top_genes <- ma_df[
  order(ma_df$padj),
][1:15, ]

#################################################
# Plot
#################################################

p_ma <- ggplot(
  ma_df,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange
  )
) +

  geom_point(
    aes(color = Status),
    alpha = 0.75,
    size = 1.8
  ) +

  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = 0,
    color = "black",
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
    title = "Meropenem-treated biofilm vs Control biofilm",
    x = expression(Log[10]~"(Mean normalized counts + 1)"),
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
  "figures/GSE167137/MA_plot_GSE167137.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 300
)
