############################################################
# Robust Rank Aggregation (RRA) Meta-analysis
# Project:
# Conserved Antibiotic Tolerance Genes in
# Pseudomonas aeruginosa Biofilms
############################################################
# Clear the workspace
rm(list = ls())

library(dplyr)
library(readr)
library(RobustRankAggreg)

############################################################
# Load significant DEGs from each study
############################################################

g1 <- read_csv(
  "results/GSE10030/GSE10030_meta_input.csv"
)

g2 <- read_csv(
  "results/GSE167137/GSE167137_meta_input.csv"
)

g3 <- read_csv(
  "results/GSE65870/GSE65870_meta_input_Ciprofloxacin_vs_Control.csv"
)

g4 <- read_csv(
  "results/GSE65870/GSE65870_meta_input_Tobramycin_vs_Control.csv"
)



g1 <- g1 %>% filter(!is.na(Gene))
g2 <- g2 %>% filter(!is.na(Gene))
g3 <- g3 %>% filter(!is.na(Gene))
g4 <- g4 %>% filter(!is.na(Gene))

############################################################
# Create ranking score
#
# Higher positive score:
# Strongly upregulated genes
#
# More negative score:
# Strongly downregulated genes
############################################################

create_rank_score <- function(df){

  df %>%
    mutate(

      adj.P.Val = ifelse(
        adj.P.Val < 1e-300,
        1e-300,
        adj.P.Val
      ),

      RankScore = logFC * (-log10(adj.P.Val))

    )

}

g1 <- create_rank_score(g1)
g2 <- create_rank_score(g2)
g3 <- create_rank_score(g3)
g4 <- create_rank_score(g4)

############################################################
# Generate ranked gene lists
############################################################

rank_up <- function(df){

  df %>%
    filter(logFC > 0) %>%
    arrange(desc(RankScore)) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    pull(Gene)

}

rank_down <- function(df){

  df %>%
    filter(logFC < 0) %>%
    arrange(RankScore) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    pull(Gene)

}

############################################################
# Build ranked lists
############################################################

up_lists <- list(

  GSE10030       = rank_up(g1),
  GSE167137      = rank_up(g2),
  GSE65870_Cipro = rank_up(g3),
  GSE65870_Tobra = rank_up(g4)

)

down_lists <- list(

  GSE10030       = rank_down(g1),
  GSE167137      = rank_down(g2),
  GSE65870_Cipro = rank_down(g3),
  GSE65870_Tobra = rank_down(g4)

)

############################################################
# Remove invalid gene identifiers
############################################################

clean_rank_list <- function(x){

  x <- trimws(x)

  x <- x[!is.na(x)]

  x <- x[x != ""]

  x <- unique(x)

  return(x)

}

up_lists <- lapply(up_lists, clean_rank_list)

down_lists <- lapply(down_lists, clean_rank_list)

############################################################
# Define gene universe
############################################################

all_genes <- unique(

  unlist(
    c(
      up_lists,
      down_lists
    )
  )

)

Ngenes <- length(all_genes)

cat("Gene universe size:", Ngenes, "\n")

############################################################
# Construct rank matrices
############################################################

rmat_up <- rankMatrix(
  up_lists,
  N = Ngenes,
  full = TRUE
)

rmat_down <- rankMatrix(
  down_lists,
  N = Ngenes,
  full = TRUE
)

############################################################
# Run Robust Rank Aggregation
############################################################

rra_up <- aggregateRanks(
  rmat = rmat_up,
  method = "RRA"
)

rra_down <- aggregateRanks(
  rmat = rmat_down,
  method = "RRA"
)

############################################################
# Inspect results
############################################################

head(rra_up)

head(rra_down)

############################################################
# Select significant meta-genes
#
# Lower score indicates stronger consistency
# across studies
############################################################

core_up <- rra_up %>%
  filter(Score < 0.01)

core_down <- rra_down %>%
  filter(Score < 0.01)

sum(rra_up$Score < 0.01)

sum(rra_down$Score < 0.01)
############################################################
# Remove genes with inconsistent regulation direction
#
# These genes appear in both the upregulated and
# downregulated RRA results, indicating inconsistent
# behavior across datasets.
#
# Since the objective is to identify a conserved
# transcriptional response, genes with conflicting
# regulation patterns are excluded from downstream
# analyses.
############################################################

conflict_genes <- intersect(
  core_up$Name,
  core_down$Name
)

cat(
  "Number of conflicting genes:",
  length(conflict_genes),
  "\n"
)

print(conflict_genes)

############################################################
# Remove conflicting genes from both gene sets
############################################################

core_up <- core_up %>%
  filter(
    !Name %in% conflict_genes
  )

core_down <- core_down %>%
  filter(
    !Name %in% conflict_genes
  )

############################################################
# Merge core upregulated and downregulated genes
############################################################

core_genes <- bind_rows(

  core_up %>%
    mutate(Direction = "Up"),

  core_down %>%
    mutate(Direction = "Down")

)

############################################################
# Verify that no duplicated genes remain
############################################################

cat(
  "Remaining duplicated genes:",
  sum(duplicated(core_genes$Name)),
  "\n"
)

############################################################
# Generate unique gene list for enrichment analyses
############################################################

genes_for_enrichment <- unique(
  core_genes$Name
)

cat(
  "Final number of core genes:",
  length(genes_for_enrichment),
  "\n"
)


########################################################
#Filter core genes
########################################################

gene_presence <- data.frame(
  Gene = core_genes$Name
)

all_lists <- c(
  up_lists,
  down_lists
)


gene_presence$Studies <- sapply(
  gene_presence$Gene,
  function(g){

    sum(
      sapply(
        all_lists,
        function(x) g %in% x
      )
    )

  }
)

core_genes <- core_genes %>%
  left_join(
    gene_presence,
    by = c("Name" = "Gene")
  ) %>%
  filter(
    Studies >= 3
  )

head(core_genes, 50)

############################################################
# Save outputs
############################################################

dir.create(
  "results/RRA",
  showWarnings = FALSE,
  recursive = TRUE
)

write_csv(
  rra_up,
  "results/RRA/RRA_upregulated_genes.csv"
)

write_csv(
  rra_down,
  "results/RRA/RRA_downregulated_genes.csv"
)



############################################################
# Save final gene list
############################################################

write_csv(
  core_genes,
  "results/RRA/Core_RRA_genes.csv"
)

write.table(
  genes_for_enrichment,
  "results/RRA/core_genes_for_enrichment.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



############################################################
# Summary
############################################################

cat(
  "Core upregulated genes:",
  nrow(core_up),
  "\n"
)

cat(
  "Core downregulated genes:",
  nrow(core_down),
  "\n"
)

cat(
  "Total core genes:",
  nrow(core_genes),
  "\n"
)



########################################
#Statistics
########################################

head(core_genes, 20)

gene_presence %>%
  filter(Studies == 3)

table(gene_presence$Studies)


nrow(core_up)

nrow(core_down)

nrow(core_genes)

head(rra_up)

head(rra_down)

length(genes_for_enrichment)
head(genes_for_enrichment, 20)


############################################################
# Genes unique to GSE167137
############################################################

other_genes <- unique(
  c(
    g1$Gene,
    g3$Gene,
    g4$Gene
  )
)

g2_unique <- setdiff(
  unique(g2$Gene),
  other_genes
)

length(g2_unique)

head(g2_unique, 50)




############################################################
# Genes missing from GSE167137
############################################################

other_genes <- unique(
  c(
    g1$Gene,
    g3$Gene,
    g4$Gene
  )
)

missing_in_g2 <- setdiff(
  other_genes,
  unique(g2$Gene)
)

length(missing_in_g2)

head(missing_in_g2, 50)

############################################################
# Prepare gene list for STRING enrichment analysis
############################################################

genes_for_string <- unique(
  core_genes$Name
)

write.table(
  genes_for_string,
  "results/core_genes_STRING.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat(
  "Number of genes exported:",
  length(genes_for_string),
  "\n"
)

############################################################

## Figures

###########################################################

############################################################
# UpSet plot showing gene occurrence across datasets
############################################################

library(dplyr)
library(tidyr)
library(UpSetR)

############################################################
# Create gene universe
############################################################

all_genes <- unique(
  c(
    g1$Gene,
    g2$Gene,
    g3$Gene,
    g4$Gene
  )
)

############################################################
# Build presence/absence matrix
############################################################

upset_df <- data.frame(
  Gene = all_genes
)

upset_df$GSE10030 <- ifelse(
  all_genes %in% g1$Gene,
  1,
  0
)

upset_df$GSE167137 <- ifelse(
  all_genes %in% g2$Gene,
  1,
  0
)

upset_df$GSE65870_Cipro <- ifelse(
  all_genes %in% g3$Gene,
  1,
  0
)

upset_df$GSE65870_Tobra <- ifelse(
  all_genes %in% g4$Gene,
  1,
  0
)

############################################################
# Save publication-quality UpSet plot
############################################################
dir.create(
  "figures/meta_analysis",
  showWarnings = FALSE,
  recursive = TRUE
)

svg(
  "figures/meta_analysis/Figure1_UpSet_CoreGenes.svg",
  width = 10.67,
  height = 7.33
)

upset(
  upset_df,
  sets = c(
    "GSE10030",
    "GSE167137",
    "GSE65870_Cipro",
    "GSE65870_Tobra"
  ),
  keep.order = TRUE,
  order.by = "freq",
  mb.ratio = c(0.65, 0.35),
  text.scale = 1.6,
  point.size = 4,
  line.size = 1.3
)

dev.off()




############################################################
# Top 20 Upregulated Core Genes Identified by RRA
############################################################

library(ggplot2)
library(dplyr)

############################################################
# Select top upregulated genes
############################################################

top20_up <- rra_up %>%
  arrange(Score) %>%
  head(20) %>%
  mutate(
    MinusLog10Score = -log10(Score)
  )

############################################################
# Generate publication-quality bar plot
############################################################

p_up <- ggplot(
  top20_up,
  aes(
    x = reorder(Name, MinusLog10Score),
    y = MinusLog10Score
  )
) +
  geom_col(
    fill = "#B2182B",
    width = 0.8
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 20 Upregulated Core Genes",
    x = "",
    y = expression(-log[10](RRA~Score))
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/meta_analysis/Figure2A_Top20_Upregulated_RRA.svg",
  p_up,
  width = 8,
  height = 6,
  dpi = 300
)




############################################################
# Top 20 Downregulated Core Genes Identified by RRA
############################################################

library(ggplot2)
library(dplyr)

############################################################
# Select top downregulated genes
############################################################

top20_down <- rra_down %>%
  arrange(Score) %>%
  head(20) %>%
  mutate(
    MinusLog10Score = -log10(Score)
  )

############################################################
# Generate publication-quality bar plot
############################################################

p_down <- ggplot(
  top20_down,
  aes(
    x = reorder(Name, MinusLog10Score),
    y = MinusLog10Score
  )
) +
  geom_col(
    fill = "#2166AC",
    width = 0.8
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 20 Downregulated Core Genes",
    x = "",
    y = expression(-log[10](RRA~Score))
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/meta_analysis/Figure2B_Top20_Downregulated_RRA.svg",
  p_down,
  width = 8,
  height = 6,
  dpi = 300
)


