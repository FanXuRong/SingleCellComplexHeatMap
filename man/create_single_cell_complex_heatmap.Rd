% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_single_cell_complex_heatmap.R
\name{create_single_cell_complex_heatmap}
\alias{create_single_cell_complex_heatmap}
\title{Create Complex Heatmap for Single Cell Expression Data}
\usage{
create_single_cell_complex_heatmap(
  seurat_object,
  features,
  gene_classification = NULL,
  group_by = "seurat_clusters",
  idents = NULL,
  time_points_order = NULL,
  cell_types_order = NULL,
  color_range = c(-1, 0, 2),
  color_palette = NULL,
  max_circle_size = 2,
  row_fontsize = 8,
  col_fontsize = 9,
  col_name_rotation = 90,
  row_title_fontsize = 10,
  col_title_fontsize = 10,
  show_heatmap_legend = TRUE,
  show_percentage_legend = TRUE,
  legend_side = "right",
  cell_border_color = "grey80",
  split_pattern = "_",
  gene_color_palette = "Set1",
  time_color_palette = "Accent",
  celltype_color_palette = "Dark2",
  show_gene_grouping = NULL,
  show_time_annotation = TRUE,
  show_celltype_annotation = TRUE,
  split_by = "time",
  merge_legends = TRUE,
  percentage_legend_title = "Expression \%",
  percentage_legend_labels = c("0\%", "25\%", "50\%", "75\%", "100\%"),
  percentage_breaks = NULL,
  return_data = FALSE,
  save_plot = NULL,
  plot_width = 10,
  plot_height = 8,
  plot_dpi = 300,
  assay = NULL,
  slot = "scale.data",
  cluster_cells = TRUE,
  cluster_features = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_cols = "complete",
  color_palette_main = c("blue", "white", "red"),
  annotation_colors = NULL,
  show_feature_names = TRUE,
  feature_names_gp = NULL,
  legend_title = "Expression",
  gene_group_title = "Gene Group",
  time_point_title = "Time Point",
  cell_type_title = "Cell Type",
  show_cell_borders = TRUE,
  show_column_annotation = TRUE,
  gene_name_mapping = NULL,
  ...
)
}
\arguments{
\item{seurat_object}{A Seurat object containing single cell data}

\item{features}{Character vector of gene names to plot}

\item{gene_classification}{Named list where names are group labels and values are character vectors of gene names (default: NULL for no gene grouping)}

\item{group_by}{Character string specifying the metadata column to group by (default: "seurat_clusters")}

\item{idents}{Numeric or character vector specifying which cell groups to include (default: NULL for all)}

\item{time_points_order}{Character vector specifying order of time points. Only affects display order, not data filtering (default: NULL for automatic)}

\item{cell_types_order}{Character vector specifying order of cell types. Only affects display order, not data filtering (default: NULL for automatic)}

\item{color_range}{Numeric vector specifying color mapping break points for expression values. Its length must match color_palette if color_palette is a vector. (default: c(-1, 0, 2))}

\item{color_palette}{Character vector specifying colors for expression heatmap. Its length must match color_range. If NULL, a default palette (viridis or color_palette_main) is generated to match color_range length (default: NULL)}

\item{max_circle_size}{Numeric specifying maximum circle radius in mm. This applies to the highest percentage value in percentage_breaks (default: 2)}

\item{row_fontsize}{Numeric specifying row name font size (default: 8)}

\item{col_fontsize}{Numeric specifying column name font size (default: 9)}

\item{col_name_rotation}{Numeric specifying column name rotation angle (default: 90)}

\item{row_title_fontsize}{Numeric specifying row title font size (default: 10)}

\item{col_title_fontsize}{Numeric specifying column title font size (default: 10)}

\item{show_heatmap_legend}{Logical indicating whether to show heatmap legend (default: TRUE)}

\item{show_percentage_legend}{Logical indicating whether to show percentage legend (default: TRUE)}

\item{legend_side}{Character string specifying legend position (default: "right")}

\item{cell_border_color}{Character string specifying cell border color (default: "grey80")}

\item{split_pattern}{Character string used to split column names for parsing (default: "_")}

\item{gene_color_palette}{Character string specifying palette name OR character vector of colors for gene groups (default: "Set1")}

\item{time_color_palette}{Character string specifying palette name OR character vector of colors for time points (default: "Accent")}

\item{celltype_color_palette}{Character string specifying palette name OR character vector of colors for cell types (default: "Dark2")}

\item{show_gene_grouping}{Logical indicating whether to show gene grouping (default: TRUE if gene_classification provided)}

\item{show_time_annotation}{Logical indicating whether to show time point annotation (default: TRUE)}

\item{show_celltype_annotation}{Logical indicating whether to show cell type annotation (default: TRUE)}

\item{split_by}{Character string specifying how to split columns: "time", "celltype", or "none" (default: "time")}

\item{merge_legends}{Logical indicating whether to merge legends (default: TRUE)}

\item{percentage_legend_title}{Character string for percentage legend title (default: "Expression \%")}

\item{percentage_legend_labels}{Character vector for percentage legend labels}

\item{percentage_breaks}{Numeric vector specifying actual percentage values corresponding to labels}

\item{return_data}{Logical; if TRUE, return underlying data instead of drawing only}

\item{save_plot}{File path to save the drawn heatmap (PNG)}

\item{plot_width}{Numeric; width in inches for saving}

\item{plot_height}{Numeric; height in inches for saving}

\item{plot_dpi}{Numeric; resolution (DPI) for saved plot}

\item{assay}{Seurat assay name to extract data from}

\item{slot}{Seurat slot name within assay (e.g., "scale.data", "data")}

\item{cluster_cells}{Logical; whether to cluster columns (cells)}

\item{cluster_features}{Logical; whether to cluster rows (features)}

\item{clustering_distance_rows}{Distance metric for row clustering}

\item{clustering_distance_cols}{Distance metric for column clustering}

\item{clustering_method_rows}{Clustering method for rows}

\item{clustering_method_cols}{Clustering method for columns}

\item{color_palette_main}{Fallback color palette when viridis unavailable}

\item{annotation_colors}{Named list of custom annotation colors}

\item{show_feature_names}{Logical; whether to show feature (row) names}

\item{feature_names_gp}{gpar object controlling feature name appearance}

\item{legend_title}{Character; title for main heatmap legend}

\item{gene_group_title}{Character string for gene group annotation title (default: "Gene Group")}

\item{time_point_title}{Character string for time point annotation title (default: "Time Point")}

\item{cell_type_title}{Character string for cell type annotation title (default: "Cell Type")}

\item{show_cell_borders}{Logical indicating whether to show cell border lines (default: TRUE)}

\item{show_column_annotation}{Logical indicating whether to show column annotations (default: TRUE)}

\item{gene_name_mapping}{Named character vector for mapping gene names, where names are original gene names and values are display names (default: NULL)}

\item{...}{Additional arguments passed to ComplexHeatmap::Heatmap()}
}
\value{
A ComplexHeatmap object. If return_data is TRUE, returns a list containing the heatmap object and underlying data matrices.
}
\description{
Creates a complex heatmap that displays both gene expression levels (as color intensity) 
and expression percentages (as circle sizes) for single cell RNA-seq data.
This function provides extensive customization options while maintaining ease of use.
}
