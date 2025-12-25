# --------------------------
# Load Required Libraries
# --------------------------
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(ggrepel)
library(ggprism)
library(patchwork)

# --------------------------
# 1. Data Preprocessing (Simplified)
# --------------------------
# Step 1: Get immune suppression genes (core gene list)
# Note: Replace with your gene list generation logic
immune_supp_genes <- unique(c(
  # Gene selection from multiple datasets (top 10% genes)
  read.delim(file_list[1:3], sep="\t")[,1][1:ceiling(nrow(read.delim(file_list[1]))*0.1)],
  read.delim(file_list[4:10], sep="\t")[,1][1:ceiling(nrow(read.delim(file_list[4]))*0.1)]
))

# Step 2: Filter LR pairs and get final gene list
lr_data <- read.table("LR_9419.txt", sep="\t", header=TRUE)
lr_filtered_1 <- lr_data[lr_data[,1] %in% immune_supp_genes,]
lr_filtered_2 <- lr_data[lr_data[,2] %in% immune_supp_genes,]
final_gene_list <- unique(c(lr_filtered_1[,2], lr_filtered_2[,1]))

# Step 3: TCGA data statistics (simplified)
tcga_stats <- data.frame(
  cancertype = character(),
  tumor = numeric(),
  normal = numeric(),
  paired = numeric(),
  stringsAsFactors = FALSE
)

# --------------------------
# 2. Core Visualization Functions (Reusable)
# --------------------------
# Boxplot theme for consistent styling
plot_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 10, colour='black', angle = 30, hjust = 0.7),
    axis.text.y = element_text(size = 12, colour='black'),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

# Facet theme for paired sample plots
facet_theme <- theme(
  strip.background = element_blank(),
  panel.spacing.x = unit(0, "cm"),
  panel.spacing.y = unit(0, "cm"),
  strip.text.x = element_text(size=7),
  strip.text.y = element_text(size=7, face="bold"),
  axis.text.x = element_text(angle=45, vjust=1, hjust=1),
  legend.position = 'none'
)

# --------------------------
# 3. TCGA Expression Profile Visualization (All Samples)
# --------------------------
# Prepare expression data (log2 transformed)
# exp_data structure: Expression | group (Tumor/Normal) | Type (Cancer Type)
plot_tcga_all <- function(exp_data) {
  # Filter cancer types with normal samples > 4
  cancer_subset <- exp_data[exp_data$Type %in% c(
    'BLCA','BRCA','COAD','KICH','KIRP','LUAD','LUSC','PRAD','READ','UCEC',
    'CHOL','ESCA','HNSC','KIRC','LIHC','STAD'
  ),]
  
  # Main boxplot with statistical comparison
  p <- ggplot(cancer_subset, aes(
    x = Type, 
    y = Expression, 
    color = factor(group, levels = c("Tumor", "Normal")),
    fill = factor(group, levels = c("Tumor", "Normal"))
  )) +
    # Error bar
    stat_boxplot(
      geom = "errorbar",
      size = 0.6,
      width = 0.4,
      linetype = "solid",
      position = position_dodge(.7),
      color = "black"
    ) +
    # Boxplot (hide outliers)
    geom_boxplot(
      outlier.colour = NA,
      size = 0.6,
      width = 0.4,
      linetype = "solid",
      position = position_dodge(.7),
      color = "black"
    ) +
    # Color scheme (JAMA style)
    scale_fill_jama(name = "Group") +
    # Wilcox test for significance
    stat_compare_means(
      label = "p.signif",
      color = "black",
      method = "wilcox.test",
      size = 4,
      label.y = 4.5
    ) +
    # Y-axis range
    coord_cartesian(ylim = c(1, 5)) +
    # Axis labels
    labs(x = "Cancer Type", y = "Expression") +
    plot_theme
  
  return(p)
}

# Plot 1: Low expression cancers
low_expr_cancers <- exp_data[exp_data$Type %in% c(
  'BLCA','BRCA','COAD','KICH','KIRP','LUAD','LUSC','PRAD','READ','UCEC'
),]
p_low <- plot_tcga_all(low_expr_cancers)
ggsave("tcga_low_expression_cancers.pdf", p_low, width = 6, height = 3)

# Plot 2: High expression cancers
high_expr_cancers <- exp_data[exp_data$Type %in% c(
  'CHOL','ESCA','HNSC','KIRC','LIHC','STAD'
),]
p_high <- plot_tcga_all(high_expr_cancers)
ggsave("tcga_high_expression_cancers.pdf", p_high, width = 5, height = 3)

# --------------------------
# 4. TCGA Paired Sample Visualization
# --------------------------
plot_tcga_paired <- function(paired_exp_data) {
  p <- ggplot(paired_exp_data, aes(
    x = factor(group, levels = c("Tumor", "Normal")),
    y = Expression
  )) +
    # Error bar
    stat_boxplot(
      geom = "errorbar",
      size = 0.6,
      width = 0.4,
      linetype = "solid",
      position = position_dodge(.7),
      color = "black"
    ) +
    # Boxplot with paired sample lines
    geom_boxplot(
      aes(fill = factor(group, levels = c("Tumor", "Normal"))),
      outlier.colour = NA,
      size = 0.6,
      width = 0.4,
      linetype = "solid",
      position = position_dodge(.7),
      color = "black"
    ) +
    # Facet by cancer type
    facet_wrap(~ Type, scales = 'free_y', nrow = 2) +
    # JAMA color scheme
    scale_fill_jama(name = "Group") +
    # Paired sample points and lines
    geom_point(size = 1.5, fill = 'gray', shape = 21, color = 'black') +
    geom_line(aes(group = sample), color = 'grey80', lwd = 0.5) +
    # Paired t-test
    stat_compare_means(
      label = "p.signif",
      color = "black",
      method = "t.test",
      size = 4,
      paired = TRUE
    ) +
    # Axis labels
    labs(x = "Sample Type", y = "Expression") +
    plot_theme +
    facet_theme
  
  return(p)
}

# Generate paired sample plot
p_paired <- plot_tcga_paired(paired_exp_data)
ggsave("tcga_paired_samples.pdf", p_paired, width = 14, height = 3)

# --------------------------
# 5. FANTOM/GTEx Expression Visualization
# --------------------------
# Convert gene symbols to ENSEMBL IDs
gene_ids <- mapIds(org.Hs.eg.db, final_gene_list, "ENSEMBL", "SYMBOL")
gene_ids <- na.omit(gene_ids)

# FANTOM data visualization
fantom_data <- read.table("FANTOM_mRNA_cpm_34.txt")
fantom_data <- log(fantom_data + 1)
fantom_data <- fantom_data[rownames(fantom_data) %in% gene_ids,]
fantom_long <- pivot_longer(
  cbind(gene = rownames(fantom_data), fantom_data),
  cols = -gene,
  names_to = "tissue",
  values_to = "expression"
)

p_fantom <- ggplot(fantom_long, aes(x = tissue, y = expression, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_classic() +
  stat_compare_means(method = "anova", label.y = 10) +
  stat_compare_means(
    label = 'p.signif',
    method = "t.test",
    ref.group = ".all.",
    label.y = 9,
    hide.ns = TRUE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 12)
  )
ggsave("fantom_tissue_expression.pdf", p_fantom, width = 8, height = 6)

# GTEx data visualization (similar structure to FANTOM)
gtex_data <- read.table("GTEx_GeneExp.txt", header = TRUE, row.names = 1)
gtex_data <- log(gtex_data + 1)
gtex_data <- gtex_data[rownames(gtex_data) %in% gene_ids,]
gtex_long <- pivot_longer(
  cbind(gene = rownames(gtex_data), gtex_data),
  cols = -gene,
  names_to = "tissue",
  values_to = "expression"
)

p_gtex <- ggplot(gtex_long, aes(x = tissue, y = expression, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_classic() +
  stat_compare_means(method = "anova", label.y = 7) +
  stat_compare_means(
    label = 'p.signif',
    method = "t.test",
    ref.group = ".all.",
    label.y = 6.5,
    hide.ns = TRUE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 12)
  )
ggsave("gtex_tissue_expression.pdf", p_gtex, width = 8, height = 6)

# --------------------------
# 6. Immune Checkpoint Correlation Visualization
# --------------------------
# Combine GSE datasets (PD-1/PD-L1/CTLA4 correlation)
gse_combined <- rbind(
  read.table("result_total_PDCD1.txt", header = TRUE),
  read.table("result_total_PDL1.txt", header = TRUE),
  read.table("result_total_CTLA4.txt", header = TRUE)
)
colnames(gse_combined) <- c('Estimate','Pvalue','Group','Gene')
gse_combined$Estimate <- abs(gse_combined$Estimate)
gse_combined$Gene <- paste(gse_combined$Gene, str_extract(gse_combined$Gene, "GSE\\d+"))

# Set gene order for plotting
gse_combined$Gene <- factor(gse_combined$Gene, levels = c(
  'CTLA-4(GSE91061)','PD-L1(GSE91061)','PD-1(GSE91061)',
  'CTLA-4(GSE165252)','PD-L1(GSE165252)','PD-1(GSE165252)',
  'CTLA-4(GSE115821)','PD-L1(GSE115821)','PD-1(GSE115821)',
  'CTLA-4(GSE176307)','PD-L1(GSE176307)','PD-1(GSE176307)'
))

# Boxplot for immune checkpoint correlation
p_checkpoint_box <- ggplot(gse_combined, aes(x = Gene, y = Estimate, color = Group, fill = Group)) +
  stat_boxplot(
    geom = "errorbar",
    size = 0.6,
    width = 0.6,
    linetype = "solid",
    position = position_dodge(.7),
    color = "black"
  ) +
  geom_boxplot(
    outlier.colour = NA,
    size = 0.6,
    width = 0.6,
    linetype = "solid",
    position = position_dodge(.7),
    color = "black"
  ) +
  coord_flip() +
  scale_fill_jama(name = "Group") +
  stat_compare_means(
    label = "p.signif",
    color = "black",
    method = "wilcox.test",
    size = 6,
    label.y = 1
  ) +
  labs(x = NULL, y = "Correlation Coefficient (R)") +
  plot_theme
ggsave("immune_checkpoint_boxplot.pdf", p_checkpoint_box, width = 9, height = 12)

# Violin plot for immune checkpoint correlation
p_checkpoint_violin <- ggplot(gse_combined, aes(x = Gene, y = Estimate, fill = Group)) +
  geom_violin(
    color = "white",
    alpha = 0.7,
    trim = FALSE,
    scale = 'width'
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#486b98","#b93735")) +
  labs(x = NULL, title = NULL) +
  coord_flip() +
  theme(
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black')
  ) +
  stat_compare_means(label = "p.signif", label.y = 1, size = 6)
ggsave("immune_checkpoint_violin.pdf", p_checkpoint_violin, width = 6, height = 9)