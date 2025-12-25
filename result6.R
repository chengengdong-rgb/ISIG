# --------------------------
# Load Required Libraries (Visualization Focus)
# --------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(ggpubr)
library(scales)

# --------------------------
# Core Visualization Setup (Reusable Elements)
# --------------------------
# Color palettes (centralized for consistency)
correlation_colors <- scale_fill_gradient2(
  low = 'blue', 
  high = 'red',
  mid = 'white', 
  limit = c(-1, 1),
  name = paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")
)

scatter_plot_theme <- theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid = element_blank()
  )

heatmap_theme <- theme(
  axis.text.x = element_text(size = 8, angle = 30, hjust = 1, color = "black"),
  axis.text.y = element_text(size = 8, color = "black"),
  axis.ticks.y = element_blank(),
  panel.background = element_blank(),
  legend.key.size = unit(1, 'cm'),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10)
)

survival_plot_palette <- c("#E7B800", "#2E9FDF")

# --------------------------
# 1. Core Visualization Functions (Reusable)
# --------------------------
# Plot correlation scatter plot (Pearson/Spearman)
plot_gene_correlation <- function(data, x_var, y_var, method = "pearson", title = "", output_file) {
  # Set correlation coefficient label
  cor_label <- ifelse(method == "pearson", "R", "rho")
  
  p <- ggplot(data, aes(!!sym(x_var), !!sym(y_var))) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
    scatter_plot_theme +
    xlab(title) +
    ylab('log(TMB)') +
    stat_cor(
      method = method, 
      aes(x = !!sym(x_var), y = !!sym(y_var)),
      cor.coef.name = cor_label,
      digits = 2
    )
  
  ggsave(output_file, p, width = 7, height = 6)
  return(p)
}

# Plot correlation heatmap with significance stars
plot_correlation_heatmap <- function(cor_data, output_file) {
  # Prepare significance labels
  cor_data <- cor_data %>%
    mutate(
      star = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ ""
      ),
      p = ifelse(is.na(p), 1, p)
    )
  
  p <- ggplot(cor_data, aes(cancerTYpe, cellType)) +
    geom_tile(aes(fill = R)) +
    geom_text(aes(label = star), color = "black", size = 4) +
    correlation_colors +
    labs(x = NULL, y = NULL) +
    heatmap_theme
  
  ggsave(output_file, p, width = 12, height = 8)
  return(p)
}

# Plot survival analysis (KM curve)
plot_survival_km <- function(data, time_var, event_var, group_var, title = "", output_file) {
  # Create survival object
  surv_obj <- Surv(data[[time_var]], data[[event_var]])
  
  # Calculate optimal cutpoint
  cutpoint <- surv_cutpoint(
    data, 
    time = time_var, 
    event = event_var, 
    variables = group_var,
    minprop = 0.1
  )
  
  # Add high/low expression group
  data <- data %>%
    mutate(
      gene_expression = ifelse(!!sym(group_var) > cutpoint[["cutpoint"]][1, 1], "high", "low")
    )
  
  # Fit KM model
  km_fit <- survfit(surv_obj ~ gene_expression, data = data)
  
  # Generate plot
  p <- ggsurvplot(
    km_fit,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    title = title,
    palette = survival_plot_palette,
    ggtheme = theme_bw()
  )
  
  ggsave(output_file, p$plot, width = 8, height = 6)
  return(p)
}

# --------------------------
# 2. Data Aggregation Function (Minimal)
# --------------------------
# Aggregate correlation results across cancer types
aggregate_correlation_results <- function(file_paths, cancer_types, value_col = 3) {
  # Initialize results dataframe
  combined_results <- read.table(file_paths[1], header = TRUE, sep = "\t") %>%
    select(all_of(value_col))
  
  colnames(combined_results) <- cancer_types[1]
  rownames(combined_results) <- rownames(read.table(file_paths[1], header = TRUE, sep = "\t"))
  
  # Combine all cancer types
  for (i in 2:length(file_paths)) {
    temp_data <- read.table(file_paths[i], header = TRUE, sep = "\t") %>%
      select(all_of(value_col))
    
    colnames(temp_data) <- cancer_types[i]
    combined_results <- cbind(combined_results, temp_data)
  }
  
  return(combined_results)
}

# --------------------------
# 3. Main Visualization Workflow (Minimal Data Prep)
# --------------------------
# --------------------------
# 3.1 Gene vs TMB Correlation Plots
# --------------------------
# Example: CD53 vs TMB (replace with your preprocessed data)
# Assume preprocessed data is loaded as a list by cancer type
# cd53_tmb_data <- readRDS("cd53_tmb_combined_data.rds")

# Generate correlation plots for all cancer types
# for (cancer_type in names(cd53_tmb_data)) {
#   # Pearson correlation
#   plot_gene_correlation(
#     data = cd53_tmb_data[[cancer_type]],
#     x_var = "gene",
#     y_var = "total_perMB_log",
#     method = "pearson",
#     title = cancer_type,
#     output_file = paste0("CD53_", cancer_type, "_pearson.pdf")
#   )
#   
#   # Spearman correlation
#   plot_gene_correlation(
#     data = cd53_tmb_data[[cancer_type]],
#     x_var = "gene",
#     y_var = "total_perMB_log",
#     method = "spearman",
#     title = cancer_type,
#     output_file = paste0("CD53_", cancer_type, "_spearman.pdf")
#   )
# }

# --------------------------
# 3.2 Immune Infiltration Correlation Heatmap
# --------------------------
# Load precomputed correlation results (minimal data prep)
# correlation_r <- aggregate_correlation_results(
#   file_paths = list.files("ITGA3_correlation_results/", pattern = "*_r.txt", full.names = TRUE),
#   cancer_types = gsub("_r.txt", "", list.files("ITGA3_correlation_results/", pattern = "*_r.txt")),
#   value_col = 3  # Spearman R column
# )
# 
# correlation_p <- aggregate_correlation_results(
#   file_paths = list.files("ITGA3_correlation_results/", pattern = "*_p.txt", full.names = TRUE),
#   cancer_types = gsub("_p.txt", "", list.files("ITGA3_correlation_results/", pattern = "*_p.txt")),
#   value_col = 4  # Spearman P column
# )

# Reshape data for heatmap
# correlation_data <- correlation_r %>%
#   mutate(cellType = rownames(.)) %>%
#   gather(key = "cancerTYpe", value = "R", -cellType) %>%
#   left_join(
#     correlation_p %>%
#       mutate(cellType = rownames(.)) %>%
#       gather(key = "cancerTYpe", value = "p", -cellType),
#     by = c("cellType", "cancerTYpe")
#   )

# Generate heatmap
# plot_correlation_heatmap(
#   cor_data = correlation_data,
#   output_file = "ITGA3_immune_infiltration_correlation_heatmap.pdf"
# )

# --------------------------
# 3.3 Survival Analysis Plots
# --------------------------
# Example: CD53 survival analysis (replace with your preprocessed data)
# survival_data <- readRDS("cd53_survival_combined_data.rds")

# Generate survival plots for all cancer types
# for (cancer_type in names(survival_data)) {
#   plot_survival_km(
#     data = survival_data[[cancer_type]],
#     time_var = "OS.time",
#     event_var = "OS",
#     group_var = "gene",
#     title = cancer_type,
#     output_file = paste0("CD53_", cancer_type, "_survival.pdf")
#   )
# }

# --------------------------
# 4. Quick Execution Example (Full Workflow)
# --------------------------
# This section demonstrates how to use the functions with sample data
# Replace with your actual data paths and parameters

# 4.1 Example: Single gene correlation plot
# sample_data <- data.frame(
#   gene = rnorm(100),
#   total_perMB_log = rnorm(100)
# )
# 
# plot_gene_correlation(
#   data = sample_data,
#   x_var = "gene",
#   y_var = "total_perMB_log",
#   method = "pearson",
#   title = "BRCA",
#   output_file = "example_correlation_plot.pdf"
# )

# 4.2 Example: Correlation heatmap
# sample_heatmap_data <- expand.grid(
#   cancerTYpe = c("ACC", "BRCA", "LUAD"),
#   cellType = c("CD8+ T cells", "CD4+ T cells", "B cells")
# ) %>%
#   mutate(
#     R = runif(n(), -1, 1),
#     p = runif(n(), 0, 1)
#   )
# 
# plot_correlation_heatmap(
#   cor_data = sample_heatmap_data,
#   output_file = "example_heatmap.pdf"
# )

# 4.3 Example: Survival plot
# sample_survival_data <- data.frame(
#   OS.time = rexp(100),
#   OS = sample(c(0,1), 100, replace = TRUE),
#   gene = rnorm(100)
# )
# 
# plot_survival_km(
#   data = sample_survival_data,
#   time_var = "OS.time",
#   event_var = "OS",
#   group_var = "gene",
#   title = "Example Cancer",
#   output_file = "example_survival_plot.pdf"
# )