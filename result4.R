# --------------------------
# Load Required Libraries (Visualization Focus)
# --------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
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
  "#FF9E4A", "#AD8BC9", "#AEC7E8", "grey", "#FF9E4A", "#ED665D", "#6DCCDA",
  "#AD8BC9", "#FFBB78", "#ED97CA", "#729ECE", "#98DF8A", "#FF9896", "#9EDAE5", 
  "#F7B6D2", "#A8786E"
)

special_cluster_rnr_colors <- c(
  "#AD8BC9", "#AEC7E8", "#FF9E4A", "grey", "#FF9E4A", "#ED665D", "#6DCCDA",
  "#AD8BC9", "#FFBB78", "#ED97CA", "#729ECE", "#98DF8A", "#FF9896", "#9EDAE5", 
  "#F7B6D2", "#A8786E"
)

response_colors <- c("#fcbbad", "#fe7aab")
dotplot_colors <- c("#ffffff", "firebrick3")

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
# 1. Core Visualization Functions (Reusable)
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
  
  return(p)
}

# Plot UMAP for special clusters (5,6 highlighted)
plot_umap_special_clusters <- function(seurat_obj, rnr_labeled = FALSE, output_file) {
  # Define color palette
  color_palette <- if (rnr_labeled) special_cluster_rnr_colors else special_cluster_colors
  
  p <- DimPlot(
    seurat_obj, 
    reduction = "umap",
    group.by = "Spcecial_cluster", 
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

# Plot cluster proportion by response (stacked percentage bar plot)
plot_cluster_proportion <- function(meta_data, output_file) {
  meta_data$seurat_clusters <- as.factor(meta_data$seurat_clusters)
  
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
  
  return(p)
}

# Plot marker gene dot plot
plot_marker_dotplot <- function(seurat_obj, marker_genes, output_file) {
  p <- DotPlot(seurat_obj, features = marker_genes) +
    RotatedAxis() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12)
    ) +
    labs(x = NULL, y = NULL) +
    scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2])
  
  pdf(output_file, height = 2, width = 8.5)
  print(p)
  dev.off()
  
  return(p)
}

# --------------------------
# 2. Combined Multi-Panel Visualization
# --------------------------
plot_combined_umap <- function(p_clusters, p_response, p_proportion, output_file) {
  p_combined <- p_clusters + p_response + p_proportion
  ggsave(output_file, p_combined, width = 18, height = 6)
  return(p_combined)
}

# --------------------------
# 3. Main Visualization Workflow (Minimal Data Prep)
# --------------------------
# Load preprocessed Seurat objects (assumes data processing is complete)
cd8_seurat <- readRDS("cd8_immune_suppression_clusters.rds")  # Immune suppression gene object
cd8_full_seurat <- readRDS("cd8_full_gene_clusters.rds")      # Full gene object for dot plots

# Load metadata (for special cluster labeling)
meta_data <- read.csv("CD8T_metadata.csv", row.names = 1)

# Add Special_cluster metadata (minimal code for visualization)
# Special clusters: 5,6 highlighted (others = "other")
meta_data$Spcecial_cluster <- meta_data$seurat_clusters
meta_data[!meta_data$Spcecial_cluster %in% c(5,6), "Spcecial_cluster"] <- "other"
cd8_seurat@meta.data$Spcecial_cluster <- meta_data$Spcecial_cluster

# --------------------------
# 4. Generate All Visualizations (One-Click Execution)
# --------------------------
# 4.1 Core UMAP Plots
p_umap_labeled <- plot_umap_clusters(
  seurat_obj = cd8_seurat,
  labeled = TRUE,
  output_file = "CD8_cluster_umap_labeled.pdf"
)

p_umap_unlabeled <- plot_umap_clusters(
  seurat_obj = cd8_seurat,
  labeled = FALSE,
  output_file = "CD8_cluster_umap_unlabeled.pdf"
)

# 4.2 Special Cluster UMAPs
p_umap_special <- plot_umap_special_clusters(
  seurat_obj = cd8_seurat,
  rnr_labeled = FALSE,
  output_file = "CD8_special_cluster_umap.pdf"
)

p_umap_special_rnr <- plot_umap_special_clusters(
  seurat_obj = cd8_seurat,
  rnr_labeled = TRUE,
  output_file = "CD8_special_cluster_RNR_umap.pdf"
)

# 4.3 Combined Multi-Panel Figure (UMAP + Response + Proportion)
p_umap_response <- plot_umap_response(cd8_seurat, "CD8_response_umap.pdf")
p_cluster_proportion <- plot_cluster_proportion(cd8_seurat@meta.data, "CD8_cluster_proportion.pdf")

p_combined <- plot_combined_umap(
  p_clusters = p_umap_unlabeled,
  p_response = p_umap_response,
  p_proportion = p_cluster_proportion,
  output_file = "CD8_combined_umap_plots.pdf"
)

# 4.4 Marker Gene Dot Plot
# Define marker genes (customize for your analysis)
marker_genes <- c(
  # Add your filtered marker genes here (replace with actual gene list)
  "CCL5", "CXCR3", "GZMK", "GZMH", "PRF1", "GNLY", "IFNG", 
  "TBX21", "EOMES", "HAVCR2", "PDCD1", "LAG3", "TIGIT"
)

# Set special cluster identities for dot plot
meta_data$Spcecial_cluster <- paste0("c", meta_data$seurat_clusters)
meta_data[!meta_data$seurat_clusters %in% c(5,6), "Spcecial_cluster"] <- "other"
Idents(cd8_full_seurat) <- meta_data$Spcecial_cluster

# Generate dot plot
p_marker_dotplot <- plot_marker_dotplot(
  seurat_obj = cd8_full_seurat,
  marker_genes = marker_genes,
  output_file = "CD8_marker_gene_dotplot.pdf"
)