## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 8,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# # Install from GitHub
# devtools::install_github("xuecheng328/SingleCellComplexHeatMap")
# 
# # Or install locally
# devtools::install_local("/path/to/SingleCellComplexHeatMap")

## ----setup--------------------------------------------------------------------
library(SingleCellComplexHeatMap)
library(Seurat)
library(dplyr)

# For this vignette, we'll use the built-in pbmc_small dataset
data("pbmc_small", package = "SeuratObject")
seurat_obj <- pbmc_small

# Add example metadata for demonstration
set.seed(123)
seurat_obj$timepoint <- sample(c("Mock", "6hpi", "24hpi"), ncol(seurat_obj), replace = TRUE)
seurat_obj$celltype <- sample(c("T_cell", "B_cell", "Monocyte"), ncol(seurat_obj), replace = TRUE)
seurat_obj$group <- paste(seurat_obj$timepoint, seurat_obj$celltype, sep = "_")

# Check the structure
head(seurat_obj@meta.data)

## ----basic_example, fig.width=8, fig.height=6---------------------------------
# Define genes to visualize
features <- c("CD3D", "CD3E", "CD8A", "IL32", "CD79A")

# Create basic heatmap
heatmap_basic <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  group_by = "group"
)

## ----advanced_example, fig.width=10, fig.height=8-----------------------------
# Define gene groups
gene_groups <- list(
  "T_cell_markers" = c("CD3D", "CD3E", "CD8A", "IL32"),
  "B_cell_markers" = c("CD79A", "CD79B", "MS4A1"),
  "Activation_markers" = c("GZMK", "CCL5")
)

# Get all genes from groups
all_genes <- c("CD3D", "CD3E", "CD8A", "IL32","CD79A", "CD79B", "MS4A1","GZMK", "CCL5")

# Create advanced heatmap with gene grouping
heatmap_advanced <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = all_genes,
  gene_classification = gene_groups,
  group_by = "group",
  time_points_order = c("Mock", "6hpi", "24hpi"),
  cell_types_order = c("T_cell", "B_cell", "Monocyte"),
  color_range = c(-2, 0, 2),
  color_palette = c("navy", "white", "firebrick"),
  max_circle_size = 3,
  split_by = "time"
)

## ----data_prep----------------------------------------------------------------
# Example of proper data preparation
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(
    # Create combined group for time course + cell type analysis
    time_celltype = paste(timepoint, celltype, sep = "_"),
    
    # Or create other combinations as needed
    cluster_time = paste(RNA_snn_res.0.8, timepoint, sep = "_")
  )

head(seurat_obj@meta.data[, c("timepoint", "celltype", "time_celltype","cluster_time")])


## ----color_customization, fig.width=10, fig.height=8--------------------------
# Custom color schemes
heatmap_colors <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = gene_groups,
  group_by = "group",
  color_range = c(-1.5, 0, 1.5,3),  # 4-point gradient
  color_palette = c("darkblue", "blue", "white", "red"), # 4 colors
  gene_color_palette = "Spectral",
  time_color_palette = "Set2",
  celltype_color_palette = "Pastel1"
)

## ----font_customization, fig.width=10, fig.height=8---------------------------
# Publication-ready styling
heatmap_publication <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = all_genes,
  gene_classification = gene_groups,
  group_by = "group",
  max_circle_size = 2.5,
  row_fontsize = 12,
  col_fontsize = 10,
  row_title_fontsize = 14,
  col_title_fontsize = 12,
  percentage_legend_title = "Fraction of cells",
  percentage_legend_labels = c("0", "20", "40", "60", "80"),
  legend_side = "right"
)

## ----clustering_control, fig.width=10, fig.height=8---------------------------
# Custom clustering
heatmap_clustering <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  group_by = "group",
  cluster_cells = TRUE,
  cluster_features = TRUE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_distance_cols = "euclidean",
  clustering_method_cols = "complete"
)

## ----time_course, fig.width=12, fig.height=8----------------------------------
# Time course focused analysis
heatmap_time <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  group_by = "group",
  time_points_order = c("Mock", "6hpi", "24hpi"),
  cell_types_order = c("T_cell", "B_cell", "Monocyte"),
  split_by = "time",
  show_celltype_annotation = TRUE,
  show_time_annotation = TRUE
)

## ----cell_type, fig.width=10, fig.height=8------------------------------------
# Cell type focused analysis
heatmap_celltype <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  group_by = "celltype",
  split_by = "celltype",
  show_time_annotation = FALSE,
  show_celltype_annotation = TRUE
)

## ----simple, fig.width=8, fig.height=6----------------------------------------
# Simple analysis
heatmap_simple <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = NULL,  # No gene grouping
  group_by = "group",
  show_time_annotation = FALSE,
  show_celltype_annotation = FALSE,
  split_by = "none"
)

## ----helper_1-----------------------------------------------------------------
# Prepare matrices
matrices <- prepare_expression_matrices(
  seurat_object = seurat_obj,
  features = features,
  group_by = "group",
  idents = NULL  # Use all groups
)

# Check the structure
dim(matrices$exp_mat)
dim(matrices$percent_mat)
head(matrices$dotplot_data)

## ----helper_2-----------------------------------------------------------------
# Create gene annotations
if (!is.null(gene_groups)) {
  gene_ann <- create_gene_annotations(
    exp_mat = matrices$exp_mat,
    percent_mat = matrices$percent_mat,
    gene_classification = gene_groups,
    color_palette = "Set1"
  )
  
  # Check results
  dim(gene_ann$exp_mat_ordered)
  levels(gene_ann$annotation_df$GeneGroup)
}

## ----helper_3-----------------------------------------------------------------
# Create cell annotations
cell_ann <- create_cell_annotations(
  exp_mat = matrices$exp_mat,
  percent_mat = matrices$percent_mat,
  time_points_order = c("Mock", "6hpi", "24hpi"),
  cell_types_order = c("T_cell", "B_cell", "Monocyte"),
  show_time_annotation = TRUE,
  show_celltype_annotation = TRUE
)

# Check results
dim(cell_ann$exp_mat_ordered)
head(cell_ann$annotation_df)

## ----save_plot, eval=FALSE----------------------------------------------------
# # Save plot to file
# heatmap_saved <- create_single_cell_complex_heatmap(
#   seurat_object = seurat_obj,
#   features = features,
#   gene_classification = gene_groups,
#   group_by = "group",
#   save_plot = "my_heatmap.png",
#   plot_width = 12,
#   plot_height = 10,
#   plot_dpi = 300
# )

## ----session_info-------------------------------------------------------------
sessionInfo()

