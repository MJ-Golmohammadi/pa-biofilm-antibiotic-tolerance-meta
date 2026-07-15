############################################################
# Figure 4: Hub Gene Network
# Publication-quality Top 50 Hub Gene Network
# Pseudomonas aeruginosa meta-analysis
############################################################

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(scales)
theme_set(theme_bw())
############################################################
# Import files
############################################################

network <- read.delim(
  "results/String Files/network.tsv",
  header = TRUE
)

hub <- read.csv(
  "results/Hubba top 50 nodes rank.csv",
  skip = 1,
  header = TRUE
)

protein <- read.delim(
  "results/String Files/protein_annotations.tsv",
  header = TRUE
)

functional <- read.delim(
  "results/String Files/functional_annotations.tsv",
  header = TRUE
)

############################################################
# Prepare hub genes
############################################################

hub <- read.csv(
  "results/Hubba top 50 nodes rank.csv",
  skip = 1,
  header = TRUE,
  check.names = FALSE
)

hub <- hub %>%
  rename(
    Gene = Name,
    MCC = Score
  )
hub_genes <- hub$Gene

############################################################
# Extract Top50 network
############################################################
network <- network %>%
  rename(
    protein1 = X.node1,
    protein2 = node2
  )

network50 <- network %>%
  filter(
    protein1 %in% hub_genes &
    protein2 %in% hub_genes
  )

############################################################
# Create functional classes
############################################################

annotation <- tibble(
  Gene = hub_genes
)

annotation$Class <- "Other"

############################################################
# Translation machinery
############################################################

translation_genes <- c(
  "tsf", "rimM", "trmD", "prs",
  "PA4463", "PA5470", "obg"
)

annotation$Class[
  grepl("^rpl|^rps|^rpm", annotation$Gene) |
  annotation$Gene %in% translation_genes
] <- "Translation"

############################################################
# Transcription
############################################################

annotation$Class[
  annotation$Gene %in%
    c(
      "rpoB",
      "rpoC",
      "rpoD"
    )
] <- "Transcription"

############################################################
# DNA repair / SOS
############################################################

annotation$Class[
  annotation$Gene %in%
    c(
      "recA",
      "lexA",
      "sulA"
    )
] <- "SOS response"

############################################################
# DNA replication
############################################################

annotation$Class[
  annotation$Gene %in%
    c(
      "dnaG"
    )
] <- "DNA replication"

############################################################
# Oxidative phosphorylation
############################################################

annotation$Class[
  grepl(
    "^nuo|^atp",
    annotation$Gene
  )
] <- "OxPhos"

############################################################
# Contractile injection system
############################################################

annotation$Class[
  annotation$Gene %in%
    c(
      "PA0614",
      "PA0618",
      "PA0620",
      "PA0621",
      "PA0622",
      "PA0624",
      "PA0627",
      "PA0628",
      "PA0629",
      "PA0633",
      "PA0634",
      "PA0635",
      "PA0636",
      "PA0637",
      "PA0638",
      "PA0641"
    )
] <- "Tailocin / R-pyocin"

############################################################
# Metabolism
############################################################

annotation$Class[
  annotation$Gene %in%
    c(
      "tpiA"
    )
] <- "Metabolism"

############################################################
# Merge MCC scores
############################################################

annotation <- annotation %>%
  left_join(
    hub,
    by = "Gene"
  )

############################################################
# Create plotting variable
############################################################

annotation$MCC_plot <- log10(annotation$MCC)

############################################################
# Build graph
############################################################

g <- graph_from_data_frame(
  d = network50,
  vertices = annotation,
  directed = FALSE
)

############################################################
# Convert to tidygraph
############################################################

tg <- as_tbl_graph(g)

############################################################
# Publication-quality colors
############################################################

my_colors <- c(

  "Translation" = "#D73027",

  "Transcription" = "#FC8D59",

  "SOS response" = "#f65cde",

  "DNA replication" = "#FEE090",

  "OxPhos" = "#1A9850",

  "Tailocin / R-pyocin" = "#4575B4",

  "Metabolism" = "#984EA3",

  "Other" = "grey70"

)

############################################################
# Plot Figure 4
############################################################

fig4 <- ggraph(
  tg,
  layout = "kk"
) +

  geom_edge_link(
    aes(
      width = combined_score,
      alpha = combined_score
    ),
    colour = "grey70"
  ) +

  scale_edge_width(
    range = c(
      0.05,
      0.7
    )
  ) +

  scale_edge_alpha(
    range = c(
      0.2,
      0.8
    )
  ) +

  geom_node_point(
    aes(
      size = MCC_plot,
      fill = Class
    ),
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +

  scale_size_continuous(
    range = c(7,24),
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill = "grey70",
        colour = "black"
      )
    )
  ) +

  guides(

    fill = guide_legend(

      override.aes = list(

        shape = 21,
        size = 14,
        colour = "black"

      )

    )

  ) +

  scale_fill_manual(
    values = my_colors
  ) +

  geom_node_text(
    aes(label = name),
    size = 4.2,
    fontface = "bold",
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.5,
    force = 10
) +

  theme_void() +

  theme(
    legend.position = "right",

    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),

    legend.title = element_text(
      face = "bold",
      size = 16
    ),

    legend.text = element_text(
      size = 14
    )
  ) +

  labs(
    title = "Top 50 Hub Gene Network",
    fill = "Functional category",
    size = "MCC score"
  )

############################################################
# Display figure
############################################################

fig4

############################################################
# Save high-quality figures
############################################################

ggsave(
  "figures/string/Figure4_HubNetwork.pdf",
  fig4,
  width = 12,
  height = 10
)


ggsave(
  "figures/string/Figure4_HubNetwork.svg",
  fig4,
  width = 12,
  height = 10
)









############################################################
# Hub gene heatmap across datasets
############################################################

library(tidyverse)
library(pheatmap)

############################################################
# Read meta-analysis datasets
############################################################

files <- c(

  GSE10030 = "results/GSE10030/GSE10030_meta_input.csv",

  GSE167137 = "results/GSE167137/GSE167137_meta_input.csv",

  GSE65870_Cipro = "results/GSE65870/GSE65870_meta_input_Ciprofloxacin_vs_Control.csv",

  GSE65870_Tobra = "results/GSE65870/GSE65870_meta_input_Tobramycin_vs_Control.csv"

)

############################################################
# Extract hub genes
############################################################

hub_genes <- hub$Gene

############################################################
# Read each dataset
############################################################
expr_list <- list()

for(i in seq_along(files)){

  tmp <- read.csv(files[i])

  tmp <- tmp %>%

    select(
      Gene,
      logFC
    ) %>%

    rename(
      !!names(files)[i] := logFC
    )

  expr_list[[i]] <- tmp
}
############################################################
# Merge all datasets
############################################################

expr <- reduce(
  expr_list,
  full_join,
  by = "Gene"
)

############################################################
# Keep only hub genes
############################################################

expr <- expr %>%

  filter(
    Gene %in% hub_genes
  )

############################################################
# Convert to matrix
############################################################

expr_mat <- expr %>%

  column_to_rownames(
    "Gene"
  ) %>%

  as.matrix()

############################################################
# Missing values
############################################################

expr_mat[is.na(expr_mat)] <- 0

############################################################
# Z-score normalization 
############################################################

expr_z <- t(
  scale(
    t(expr_mat)
  )
)

############################################################
# Create heatmap annotation
############################################################
annotation_heat <- annotation %>%
  select(
    Gene,
    Class
  )

annotation_heat <- as.data.frame(annotation_heat)

rownames(annotation_heat) <- annotation_heat$Gene

annotation_heat$Gene <- NULL

annotation_heat <- annotation_heat[
  rownames(expr_z),
  ,
  drop = FALSE
]

############################################################
# Publication colors
############################################################
ann_colors <- list(

  Class = c(

    "Translation" = "#D73027",

    "Transcription" = "#FC8D59",

    "SOS response" = "#f65cde",

    "DNA replication" = "#FEE090",

    "OxPhos" = "#1A9850",

    "Tailocin / R-pyocin" = "#4575B4",

    "Metabolism" = "#984EA3",

    "Other" = "grey70"

  )
)


############################################################
# Draw heatmap
############################################################
pheatmap(

  expr_z,

  scale = "none",

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  show_rownames = TRUE,

  show_colnames = TRUE,

  fontsize_row = 13,

  fontsize_col = 14,

  border_color = NA,

  annotation_row = annotation_heat,

  annotation_colors = ann_colors,

  color = colorRampPalette(

    c(
      "#2166AC",
      "white",
      "#B2182B"
    )

  )(100),

  main = "Hub Gene Expression Patterns"
)

############################################################
# Save figure
############################################################
graphics.off()

library(svglite)


svglite(

  "figures/string/Figure6_HubHeatmap.svg",

  width = 9,

  height = 12

)

pheatmap(

  expr_z,

  scale = "none",

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  annotation_row = annotation_heat,

  annotation_colors = ann_colors,

  fontsize_row = 13,

  fontsize_col = 14,

  border_color = NA,

  color = colorRampPalette(

    c(
      "#2166AC",
      "white",
      "#B2182B"
    )

  )(100),

  main = "Hub Gene Expression Patterns"

)

dev.off()

