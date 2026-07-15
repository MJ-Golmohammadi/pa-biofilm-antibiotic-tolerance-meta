############################################################
# Sensitivity Analysis
# Significant DEGs only
############################################################

s1 <- read.csv("results/GSE10030/GSE10030_significant.csv")
s2 <- read.csv("results/GSE167137/GSE167137_significant.csv")
s3 <- read.csv("results/GSE65870/GSE65870_significant_Ciprofloxacin_vs_Control.csv")
s4 <- read.csv("results/GSE65870/GSE65870_significant_Tobramycin_vs_Control.csv")

sig_lists <- list(
  unique(s1$Gene),
  unique(s2$Gene),
  unique(s3$Gene),
  unique(s4$Gene)
)

common_sig <- Reduce(intersect, sig_lists)

write.csv(
  data.frame(Gene = common_sig),
  "results/Meta_Analysis/Common_Significant_Genes.csv",
  row.names = FALSE
)

############################################################
# Overlap with RRA core genes
############################################################

core_rra <- read.csv(
  "results/Meta_Analysis/Core_Genes_RRA.csv"
)

validated_genes <- intersect(
  core_rra$Gene,
  common_sig
)

write.csv(
  data.frame(Gene = validated_genes),
  "results/Meta_Analysis/Validated_Core_Genes.csv",
  row.names = FALSE
)

cat("\nSensitivity analysis completed.\n")