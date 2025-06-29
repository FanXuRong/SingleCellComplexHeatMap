## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 12,
  fig.height = 8,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# # Install from GitHub
# devtools::install_github("FanXuRong/SingleCellComplexHeatMap")

## ----setup--------------------------------------------------------------------
library(SingleCellComplexHeatMap)
library(Seurat)
library(dplyr)
# Load optional color packages for testing
if (requireNamespace("ggsci", quietly = TRUE)) {
  library(ggsci)
}
if (requireNamespace("viridis", quietly = TRUE)) {
  library(viridis)
}

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

## ----basic_example, fig.width=12, fig.height=8--------------------------------
# Define genes to visualize
features <- c("CD3D", "CD3E", "CD8A", "IL32", "CD79A")

# Create basic heatmap
heatmap_basic <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  group_by = "group"
)

## ----advanced_example, fig.width=12, fig.height=8-----------------------------
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


## ----test_custom_titles, fig.width=12, fig.height=8---------------------------
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = gene_groups,
  group_by = "group",
  time_points_order = c("Mock", "6hpi", "24hpi"),
  
  # NEW: Custom annotation titles
  gene_group_title = "Gene Function",
  time_point_title = "Time Point", 
  cell_type_title = "Cell Type",
  
  split_by = "time"
)

## ----color_customization, fig.width=12, fig.height=8--------------------------
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

## ----test_ggsci_colors, fig.width=12, fig.height=8----------------------------
if (requireNamespace("ggsci", quietly = TRUE)) {
  heatmap_colors <- create_single_cell_complex_heatmap(
    seurat_object = seurat_obj,
    features = features,
    gene_classification = gene_groups,
    group_by = "group",
    time_points_order = c("Mock", "6hpi", "24hpi"),
    
    # NEW: Using ggsci color vectors
    gene_color_palette = pal_npg()(3),
    time_color_palette = pal_lancet()(3),
    celltype_color_palette = pal_jama()(4),
    
    # Custom expression heatmap colors
    color_range = c(-2, 0, 2),
    color_palette = c("#2166AC", "#F7F7F7", "#B2182B")
  )
} 

## ----test_viridis_custom_colors, fig.width=12, fig.height=8-------------------
if (requireNamespace("ggsci", quietly = TRUE)) {
  heatmap_colors <- create_single_cell_complex_heatmap(
    seurat_object = seurat_obj,
    features = features,
    gene_classification = gene_groups,
    group_by = "group",
    time_points_order = c("Mock", "6hpi", "24hpi"),
    
    # NEW: Using ggsci color vectors
    gene_color_palette = pal_npg()(3),
    time_color_palette = pal_lancet()(3),
    celltype_color_palette = pal_jama()(4),
    
    # Custom expression heatmap colors
    color_range = c(-2, 0, 2),
    color_palette = c("#2166AC", "#F7F7F7", "#B2182B")
  )
} 

## ----font_customization, fig.width=12, fig.height=8---------------------------
# Publication-ready styling
heatmap_publication <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = all_genes,
  gene_classification = gene_groups,
  group_by = "group",
  max_circle_size = 2.5,
  row_fontsize = 12,
  col_fontsize = 12,
  row_title_fontsize = 14,
  col_title_fontsize = 12,
  percentage_legend_title = "Fraction of cells",
  percentage_legend_labels = c("0", "20", "40", "60", "80"),
  legend_side = "right"
)

## ----test_visual_control, fig.width=12, fig.height=8--------------------------
heatmap_con <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = gene_groups,
  group_by = "group",
  
  # NEW: Visual control parameters
  show_cell_borders = FALSE,
  show_column_annotation = FALSE,
  
  # Other parameters for a clean plot
  split_by = "none",
  cluster_cells = TRUE
)

## ----test_gene_mapping, fig.width=12, fig.height=8----------------------------
# Create a mapping for a subset of genes
gene_mapping <- c(
  "CD3D" = "T-cell Receptor CD3d",
  "CD79A" = "B-cell Antigen Receptor CD79a",
  "GZMK" = "Granzyme K",
  "NKG7" = "Natural Killer Cell Granule Protein 7"
)

heatmap_map <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = gene_groups,
  group_by = "group",
  
  # NEW: Gene name mapping
  gene_name_mapping = gene_mapping,
  
  row_fontsize = 9
)

## ----clustering_control, fig.width=12, fig.height=8---------------------------
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

## ----cell_type, fig.width=12, fig.height=8------------------------------------
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
heatmap_sample <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = NULL,  # No gene grouping
  group_by = "group",
  show_time_annotation = FALSE,
  show_celltype_annotation = FALSE,
  split_by = "none"
)

## ----test_comprehensive, fig.width=14, fig.height=12--------------------------
create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = features,
  gene_classification = gene_groups,
  group_by = "group",
  time_points_order = c("Mock", "6hpi", "24hpi"),
  
  # --- New Features ---
  gene_group_title = "Functional Category",
  time_point_title = "Time Post-Infection",
  cell_type_title = "Cell Identity",
  show_cell_borders = TRUE,
  cell_border_color = "white",
  gene_name_mapping = c("MS4A1" = "CD20"),
  
  # --- Color Customization ---
  color_range = c(-2, 0, 2),
  color_palette = c("#0072B2", "white", "#D55E00"), # Colorblind-friendly
  gene_color_palette = "Dark2",
  time_color_palette = "Set2",
  celltype_color_palette = "Paired",
  
  # --- Layout and Font ---
  row_fontsize = 10,
  col_fontsize = 9,
  row_title_fontsize = 12,
  col_title_fontsize = 12,
  col_name_rotation = 45,
  legend_side = "right",
  merge_legends = TRUE,
  
  # --- Clustering and Splitting ---
  cluster_features = FALSE, # Rows are already grouped
  cluster_cells = FALSE,    # Columns are already grouped
  split_by = "time"
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

