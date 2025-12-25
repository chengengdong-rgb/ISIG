# --------------------------
# Load Required Libraries (Visualization Focus)
# --------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggsci)
library(ggpubr)
library(scales)

# --------------------------
# Core Visualization Setup (Reusable Elements)
# --------------------------
# Custom color palettes (centralized for consistency)
cluster_colors <- c(
  "#AEC7E8", "#67BF5C", "#CDCC5D", "#FF9E4A", "#ED665D", "#6DCCDA",
  "#AD8BC9", "#FFBB78", "#ED97CA", "#729ECE", "#98DF8A", "#FF9896", 
  "#9EDAE5", "#F7B6D2", "#A8786E"
)

special_cluster_colors <- c(
  "#AEC7E8", "#67BF5C", "#CDCC5D", "grey", "#FF9E4A", "#ED665D", "#6DCCDA",
  "#AD8BC9", "#FFBB78", "#ED97CA", "#729ECE", "#98DF8A", "#FF9896", "#9EDAE5", 
  "#F7B6D2", "#A8786E"
)

special_cluster_rnr_colors <- c(
  '#F8B0A5', "grey", "#F279A4", "#AEC7E8", "#67BF5C", "#CDCC5D", "grey", 
  "#FF9E4A", "#ED665D", "#6DCCDA", "#AD8BC9", "#FFBB78", "#ED97CA", 
  "#729ECE", "#98DF8A", "#FF9896", "#9EDAE5", "#F7B6D2", "#A8786E"
)

response_colors <- c("#fcbbad", "#fe7aab")
marker_colors <- c("#ebe8e8", "red")

# Standard UMAP theme (applied to all UMAP plots)
umap_theme <- theme(
  legend.title = element_blank(),
  legend.key = element_rect(fill = 'white'),
  legend.text = element_text(size = 20),
  legend.key.size = unit(1, 'cm'),
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(size = 14, angle = 0),
  axis.text.y = element_text(size = 14, angle = 0, hjust = 1)
)

# --------------------------
# 1. UMAP Visualization Functions (Reusable)
# --------------------------
# Plot UMAP by cluster (with/without labels)
plot_umap_clusters <- function(seurat_obj, labeled = TRUE, output_file) {
  p <- DimPlot(
    seurat_obj, 
    reduction = "umap",
    group.by = "seurat_clusters", 
    pt.size = 0.5,
    label = labeled,
    cols = cluster_colors
  ) +
    umap_theme +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(x = "UMAP_1", y = "UMAP_2")
  
  ggsave(output_file, p, width = 7, height = 6)
  return(p)
}

# Plot UMAP by treatment response
plot_umap_response <- function(seurat_obj, output_file) {
  p <- DimPlot(
    seurat_obj, 
    reduction = "umap",
    group.by = "Response", 
    pt.size = 0.5,
    cols = response_colors
  ) +
    umap_theme +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(x = "UMAP_1", y = "UMAP_2")
  
  ggsave(output_file, p, width = 7, height = 6)
  return(p)
}

# Plot UMAP for special clusters (2,6,7 highlighted)
plot_umap_special_clusters <- function(seurat_obj, rnr_labeled = FALSE, output_file) {
  # Define color palette
  color_palette <- if (rnr_labeled) special_cluster_rnr_colors else special_cluster_colors
  
  p <- DimPlot(
    seurat_obj, 
    reduction = "umap",
    group.by = "Special_cluster", 
    pt.size = 0.5,
    label = TRUE,
    cols = color_palette
  ) +
    umap_theme +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(x = "UMAP_1", y = "UMAP_2")
  
  ggsave(output_file, p, width = 7, height = 6)
  return(p)
}

# --------------------------
# 2. Supplementary Visualization Functions
# --------------------------
# Plot cluster proportion by response (stacked percentage bar plot)
plot_cluster_proportion <- function(meta_data, output_file) {
  p <- ggplot(data = meta_data, aes(x = Response, fill = seurat_clusters)) +
    geom_bar(width = 0.3, position = "fill") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = cluster_colors) +
    labs(y = "Percentage") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 14, angle = 0),
      axis.text.y = element_text(size = 14, angle = 0, hjust = 1),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.key.size = unit(1, 'cm')
    )
  
  ggsave(output_file, p, width = 8, height = 6)
  return(p)
}

# Plot marker gene feature plots
plot_marker_features <- function(seurat_obj, marker_genes, output_file) {
  p <- FeaturePlot(
    seurat_obj, 
    features = marker_genes,
    cols = marker_colors
  )
  
  ggsave(output_file, p, width = 8, height = 9)
  return(p)
}

# --------------------------
# 3. Combined Visualization (Multi-Panel Figure)
# --------------------------
plot_combined_umap <- function(p_clusters, p_response, p_proportion, output_file) {
  p_combined <- p_clusters + p_response + p_proportion
  ggsave(output_file, p_combined, width = 18, height = 6)
  return(p_combined)
}

# --------------------------
# 4. Main Visualization Workflow (Minimal Data Prep)
# --------------------------
# Load preprocessed Seurat object (assumes data processing is complete)
# Replace with your Seurat object path
cd4_seurat <- readRDS("cd4_tcells_preprocessed.rds")

# Add Special_cluster metadata (minimal code for visualization)
cluster_meta <- cd4_seurat@meta.data
cluster_meta$Special_cluster <- as.numeric(as.character(cluster_meta$seurat_clusters))
cluster_meta[!cluster_meta$Special_cluster %in% c(2,6,7), "Special_cluster"] <- "other"
cd4_seurat@meta.data$Special_cluster <- cluster_meta$Special_cluster

# Add R/NR labeling to special clusters (for secondary plot)
cluster_meta$Special_cluster_RNR <- cluster_meta$Special_cluster
cluster_meta[cluster_meta$Special_cluster_RNR == 2, "Special_cluster_RNR"] <- cluster_meta[cluster_meta$Special_cluster_RNR == 2, "Response"]
cluster_meta[cluster_meta$Special_cluster_RNR == 6, "Special_cluster_RNR"] <- cluster_meta[cluster_meta$Special_cluster_RNR == 6, "Response"]
cluster_meta[cluster_meta$Special_cluster_RNR == 7, "Special_cluster_RNR"] <- cluster_meta[cluster_meta$Special_cluster_RNR == 7, "Response"]
cd4_seurat@meta.data$Special_cluster_RNR <- cluster_meta$Special_cluster_RNR

# --------------------------
# 5. Generate All Visualizations (One-Click Execution)
# --------------------------
# 5.1 Core UMAP Plots
p_umap_labeled <- plot_umap_clusters(
  seurat_obj = cd4_seurat,
  labeled = TRUE,
  output_file = "CD4_cluster_umap_labeled.pdf"
)

p_umap_unlabeled <- plot_umap_clusters(
  seurat_obj = cd4_seurat,
  labeled = FALSE,
  output_file = "CD4_cluster_umap_unlabeled.pdf"
)

p_umap_response <- plot_umap_response(
  seurat_obj = cd4_seurat,
  output_file = "CD4_response_umap.pdf"
)

# 5.2 Special Cluster UMAPs
p_umap_special <- plot_umap_special_clusters(
  seurat_obj = cd4_seurat,
  rnr_labeled = FALSE,
  output_file = "CD4_special_cluster_umap.pdf"
)

p_umap_special_rnr <- plot_umap_special_clusters(
  seurat_obj = cd4_seurat,
  rnr_labeled = TRUE,
  output_file = "CD4_special_cluster_RNR_umap.pdf"
)

# 5.3 Supplementary Plots
p_cluster_proportion <- plot_cluster_proportion(
  meta_data = cluster_meta,
  output_file = "CD4_cluster_proportion.pdf"
)

# 5.4 Combined Multi-Panel Figure
p_combined <- plot_combined_umap(
  p_clusters = p_umap_unlabeled,
  p_response = p_umap_response,
  p_proportion = p_cluster_proportion,
  output_file = "CD4_combined_umap_plots.pdf"
)

# 5.5 Marker Gene Feature Plots
# Load full-gene Seurat object (for feature plotting)
cd4_full_seurat <- readRDS("cd4_tcells_full.rds")

p_marker_features <- plot_marker_features(
  seurat_obj = cd4_full_seurat,
  marker_genes = c('CCL5','CXCR3','GZMK','GZMH','VCAM1','XCL1'),
  output_file = "CD4_marker_gene_feature_plots.pdf"
)