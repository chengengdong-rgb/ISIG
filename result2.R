# --------------------------
# Load Required Libraries
# --------------------------
library(corrplot)
library(ggplot2)
library(gcookbook)
library(gg.gap)
library(ggbreak)
library(stringr)
library(pheatmap)

# --------------------------
# 1. Cluster Correlation Heatmap
# --------------------------
# Load cluster metadata (replace with your file path)
cluster_meta <- read.csv("cluster_meta.csv", row.names = 1)
cluster_meta[is.na(cluster_meta)] <- 0 # Handle missing values

# Filter clusters with FC > 2 or FC < 0.5
cluster_meta_filtered <- cluster_meta[cluster_meta$FC > 2 | cluster_meta$FC < 0.5, ]
cluster_meta_filtered <- cluster_meta_filtered[, -15] # Remove unused column
cluster_meta_t <- as.data.frame(t(cluster_meta_filtered)) # Transpose for correlation

# Define custom color palettes for heatmaps
color_palette1 <- c(
  colorRampPalette(c("#8B008B", "white"))(500),
  colorRampPalette(c("white", "#FF4500"))(500)
)
color_palette2 <- c(
  colorRampPalette(c("#11427C", "white"))(500),
  colorRampPalette(c("white", "#C31E1F"))(500)
)

# Calculate Pearson correlation and significance
cor_matrix <- cor(cluster_meta_t, method = "pearson")
cor_test <- cor.mtest(cluster_meta_t, method = "pearson", conf.level = 0.95)

# Plot 1: Correlation heatmap (circle + significance)
pdf(file = "cluster_correlation_heatmap_1.pdf", width = 12, height = 10)
corrplot(
  cor_matrix,
  method = "circle",
  tl.col = "black", tl.cex = 1.5, tl.srt = 45, tl.pos = "lt",
  p.mat = cor_test$p, diag = TRUE, type = 'lower',
  sig.level = c(0.001, 0.01, 0.05), pch.cex = 2.2, col = color_palette1,
  insig = 'label_sig', pch.col = 'black', order = 'hclust'
)
dev.off()

# Plot 2: Correlation heatmap (alternative color palette)
pdf(file = "cluster_correlation_heatmap_2.pdf", width = 12, height = 10)
corrplot(
  cor_matrix,
  method = "circle",
  tl.col = "black", tl.cex = 1.5, tl.srt = 45, tl.pos = "lt",
  p.mat = cor_test$p, diag = TRUE, type = 'lower',
  sig.level = c(0.001, 0.01, 0.05), pch.cex = 2.2, col = color_palette2,
  insig = 'label_sig', pch.col = 'black', order = 'hclust'
)
dev.off()

# --------------------------
# 2. Cluster FC Lollipop Plot
# --------------------------
# Function to calculate cluster frequency and fold change (R/NR)
calculate_cluster_fc <- function(file_path) {
  # Read cluster data
  cluster_data <- read.csv(file_path, sep = ",", row.names = 1)
  
  # Add cluster ID with sample prefix
  sample_prefix <- strsplit(basename(file_path), '_')[[1]][1]
  cluster_data$seurat_clusters <- str_c(sample_prefix, '_C', cluster_data$seurat_clusters)
  
  # Calculate cell count per cluster
  cluster_freq <- as.data.frame(table(cluster_data$seurat_clusters))
  
  # Calculate fold change (R/NR ratio)
  for (i in 1:nrow(cluster_freq)) {
    cluster_name <- cluster_freq[i, 1]
    r_ratio <- sum(cluster_data$seurat_clusters == cluster_name & cluster_data$Response == "R") / sum(cluster_data$Response == "R")
    nr_ratio <- sum(cluster_data$seurat_clusters == cluster_name & cluster_data$Response == "NR") / sum(cluster_data$Response == "NR")
    cluster_freq[i, 3] <- r_ratio / nr_ratio
  }
  return(cluster_freq)
}

# Process all cluster files (replace with your file list)
file_list <- list.files(path = "cluster_meta/", pattern = "*.csv", full.names = TRUE)
cluster_freq_combined <- calculate_cluster_fc(file_list[1])

for (i in 2:length(file_list)) {
  cluster_freq <- calculate_cluster_fc(file_list[i])
  cluster_freq_combined <- rbind(cluster_freq_combined, cluster_freq)
}

# Rename columns and preprocess FC values
colnames(cluster_freq_combined) <- c('Cluster', 'Cell_number', 'Cluster_FC')
cluster_freq_combined$Cluster_FC <- log2(cluster_freq_combined$Cluster_FC)
cluster_freq_combined[cluster_freq_combined$Cluster_FC == '-Inf', 3] <- -2 # Handle infinite values

# Order clusters by FC (ascending)
cluster_freq_combined <- cluster_freq_combined[order(cluster_freq_combined$Cluster_FC, decreasing = FALSE), ]
cluster_freq_combined$Cluster <- factor(cluster_freq_combined$Cluster, levels = cluster_freq_combined$Cluster)

# Plot 1: Vertical lollipop plot (with broken x-axis)
p_lollipop_vertical <- ggplot(cluster_freq_combined, aes(x = Cluster_FC, y = Cluster)) +
  geom_col(width = 0.1, aes(fill = Cell_number)) + # Thin bar (stick)
  scale_x_break(c(1.5, 2.5), scales = 0.1, space = 0.1) + # Break x-axis
  geom_point(aes(size = abs(Cluster_FC), color = Cell_number)) + # Bubble (lollipop head)
  geom_vline(xintercept = c(-1, 1), lty = 8, size = 0.3, color = 'black') + # Threshold lines
  scale_size_continuous(range = c(2, 7)) + # Adjust bubble size range
  scale_fill_gradient2(low = "#486b98", mid = "#f5f2b1", high = "#b93735", midpoint = 1500) +
  scale_color_gradient2(low = "#486b98", mid = "#f5f2b1", high = "#b93735", midpoint = 1500) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

ggsave(p_lollipop_vertical, filename = "cluster_fc_lollipop_vertical.pdf", width = 7, height = 11)

# Plot 2: Horizontal lollipop plot (flipped coordinates)
p_lollipop_horizontal <- ggplot(cluster_freq_combined, aes(x = Cluster_FC, y = Cluster)) +
  geom_col(width = 0.1, aes(fill = Cell_number)) +
  scale_x_break(c(1.5, 2.5), scales = 0.1, space = 0.1) +
  geom_point(aes(size = abs(Cluster_FC), color = Cell_number)) +
  geom_vline(xintercept = c(-1, 1), lty = 8, size = 0.3, color = 'black') +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_gradient2(low = "#486b98", mid = "#f5f2b1", high = "#b93735", midpoint = 1500) +
  scale_color_gradient2(low = "#486b98", mid = "#f5f2b1", high = "#b93735", midpoint = 1500) +
  theme_classic() +
  coord_flip() + # Flip axes for horizontal plot
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

ggsave(p_lollipop_horizontal, filename = "cluster_fc_lollipop_horizontal.pdf", width = 15, height = 5)

# --------------------------
# 3. Immune Gene Expression Heatmap
# --------------------------
# Load pseudo-single cell expression data
immune_exp_data <- read.table(
  "immune_gene_pseudosinglecell_expression.txt",
  sep = "\t", header = TRUE, row.names = 1
)

# Normalize expression by cell count per cluster (replace with your normalization logic)
normalization_factors <- c(1610, 647, 211, 881, 727, 66, 147, 45, 42, 41, 39, 182, 39)
for (i in 1:length(normalization_factors)) {
  immune_exp_data[, i] <- immune_exp_data[, i] / normalization_factors[i]
}

# Define immune gene list
immune_genes <- c(
  'CCL5','CXCL13','MIF','CXCL2','CXCL16',
  'CD27','CD52','CD96','CTLA4','CXCR4','HMOX1','ICOS','KLRC1','LAG3','PDCD1','REL','RRM2','TGFB1','TIGIT','TNFRSF18','TNFRSF1B','TNFRSF4','TNFRSF9',
  'CD74','FTH1','BTG1','CD7','CD3D','CD52','HLA-DRB1','IL32','S100A4','S100A11',
  'CD69','CCR7','HLA-DPA1','HLA-DRB5','FOXP3',
  'CD28','CD47','HAVCR2','CD70','HLA-A','HLA-B','HLA-C','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRA','HLA-E'
)

# Filter immune genes and remove unwanted rows
immune_exp_filtered <- immune_exp_data[rownames(immune_exp_data) %in% immune_genes, ]
immune_exp_filtered <- immune_exp_filtered[-c(11, 16, 17, 33), ]

# Define heatmap color palette
heatmap_colors <- c(
  colorRampPalette(c("#729ECE", "white"))(50),
  colorRampPalette(c("white", "#ED665D"))(950)
)

# Plot immune gene expression heatmap
pdf("immune_gene_expression_heatmap.pdf", width = 15, height = 5)
pheatmap(
  as.data.frame(t(immune_exp_filtered)),
  cluster_row = FALSE, cluster_col = FALSE, border_color = 'gray90',
  show_rownames = TRUE, color = heatmap_colors, angle_col = '45'
)
dev.off()