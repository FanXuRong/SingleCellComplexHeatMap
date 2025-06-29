# Suppress R CMD check notes for dplyr column names
utils::globalVariables(c("TimePoint", "CellType", ".data"))

#' Create Cell Type and Time Point Annotations for Heatmap Columns
#'
#' Parses column names to extract time points and cell types, creates annotations and reorders matrices.
#'
#' @param exp_mat Expression matrix with samples as columns
#' @param percent_mat Percentage matrix with samples as columns
#' @param split_pattern Character string used to split column names (default: "_")
#' @param time_position Integer indicating position of time point in split names (default: 1)
#' @param celltype_start Integer indicating starting position of cell type in split names (default: 2)
#' @param time_points_order Character vector specifying order of time points (default: NULL for automatic)
#' @param cell_types_order Character vector specifying order of cell types (default: NULL for automatic)
#' @param time_color_palette Character string specifying palette name OR character vector of colors for time points (default: "Accent")
#' @param celltype_color_palette Character string specifying palette name OR character vector of colors for cell types (default: "Dark2")
#' @param show_time_annotation Logical indicating whether to show time point annotation (default: TRUE)
#' @param show_celltype_annotation Logical indicating whether to show cell type annotation (default: TRUE)
#' @param time_point_title Character string for time point annotation title (default: "Time Point")
#' @param cell_type_title Character string for cell type annotation title (default: "Cell Type")
#'
#' @return A list containing exp_mat_ordered (reordered expression matrix), percent_mat_ordered (reordered percentage matrix), col_annotation (ComplexHeatmap column annotation object), col_split_factor (factor for column splitting based on time points), and annotation_df (data frame with column annotations).
#'
#' @export
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils globalVariables
#'
#' @examples
#' # Load a small example Seurat object
#' data("pbmc_small", package = "SeuratObject")
#' pbmc_small$timepoint <- sample(c("0h", "6h"), ncol(pbmc_small), replace = TRUE)
#' pbmc_small$timepoint_celltype <- paste(pbmc_small$timepoint, pbmc_small$RNA_snn_res.0.8, sep = "_")
#' features <- c("CD3D", "CD79A", "MS4A1")
#' 
#' # Prepare expression matrices first
#' matrices <- prepare_expression_matrices(pbmc_small, features, group_by = "timepoint_celltype")
#' 
#' # Create cell annotations with custom ordering
#' col_annotations <- create_cell_annotations(
#'   exp_mat = matrices$exp_mat,
#'   percent_mat = matrices$percent_mat,
#'   split_pattern = "_",
#'   time_points_order = c("0h", "6h"),
#'   cell_types_order = levels(pbmc_small$RNA_snn_res.0.8)
#' )
#' 
#' # Access results
#' ordered_exp_mat <- col_annotations$exp_mat_ordered
#' 
#' @seealso \code{\link{create_single_cell_complex_heatmap}}, \code{\link{prepare_expression_matrices}}

create_cell_annotations <- function(exp_mat, 
                                  percent_mat,
                                  split_pattern = "_",
                                  time_position = 1,
                                  celltype_start = 2,
                                  time_points_order = NULL,
                                  cell_types_order = NULL,
                                  time_color_palette = "Accent",
                                  celltype_color_palette = "Dark2",
                                  show_time_annotation = TRUE,
                                  show_celltype_annotation = TRUE,
                                  time_point_title = "Time Point",
                                  cell_type_title = "Cell Type") {
  
  col_ids <- colnames(exp_mat)
  
  # Only parse if annotations are needed
  if (show_time_annotation || show_celltype_annotation) {
    parsed_parts <- strsplit(col_ids, split_pattern, fixed = TRUE)
    
    # Extract time points only if needed
    time_points_extracted <- NULL
    if (show_time_annotation) {
      time_points_extracted <- sapply(parsed_parts, `[`, time_position)
    }
    
    # Extract cell types only if needed
    cell_types_extracted <- NULL
    if (show_celltype_annotation) {
      cell_types_extracted <- sapply(parsed_parts, function(x) {
        if (length(x) >= celltype_start) {
          paste(x[celltype_start:length(x)], collapse = split_pattern)
        } else {
          # // If no underscore found, use the whole string as cell type
          x[1]
        }
      })
    }
    
    # // Set default orders if not provided
    if (show_time_annotation && is.null(time_points_order)) {
      time_points_order <- unique(time_points_extracted[!is.na(time_points_extracted)])
    }
    if (show_celltype_annotation && is.null(cell_types_order)) {
      if (is.null(cell_types_extracted)) {
        # // If no parsing needed, use column names directly as cell types
        cell_types_order <- unique(col_ids)
        cell_types_extracted <- col_ids
      } else {
        cell_types_order <- unique(cell_types_extracted[!is.na(cell_types_extracted)])
      }
    }
    
    # // Create annotation data frame with custom titles
    annotation_col_df <- data.frame(id = col_ids)
    
    if (show_time_annotation) {
      annotation_col_df[[time_point_title]] <- factor(time_points_extracted, levels = time_points_order)
    }
    
    if (show_celltype_annotation) {
      annotation_col_df[[cell_type_title]] <- factor(cell_types_extracted, levels = cell_types_order)
    }
    
    # // Arrange data frame based on what annotations are shown
    if (show_time_annotation && show_celltype_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(.data[[time_point_title]], .data[[cell_type_title]])
    } else if (show_time_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(.data[[time_point_title]])
    } else if (show_celltype_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(.data[[cell_type_title]])
    }
  } else {
    # // If no annotations are shown, just keep original order
    annotation_col_df <- data.frame(id = col_ids)
  }
  
  # // Reorder matrices
  exp_mat_ordered <- exp_mat[, annotation_col_df$id, drop = FALSE]
  percent_mat_ordered <- percent_mat[, annotation_col_df$id, drop = FALSE]
  
  # // Create color schemes with simplified approach
  annotation_list <- list()
  color_list <- list()
  
  if (show_time_annotation && time_point_title %in% colnames(annotation_col_df)) {
    n_timepoints <- length(levels(annotation_col_df[[time_point_title]]))
    time_colors <- get_colors_simplified(n_timepoints, time_color_palette)
    time_point_colors <- setNames(time_colors, levels(annotation_col_df[[time_point_title]]))
    annotation_list[[time_point_title]] <- annotation_col_df[[time_point_title]]
    color_list[[time_point_title]] <- time_point_colors
  }
  
  if (show_celltype_annotation && cell_type_title %in% colnames(annotation_col_df)) {
    n_celltypes <- length(levels(annotation_col_df[[cell_type_title]]))
    cell_colors <- get_colors_simplified(n_celltypes, celltype_color_palette)
    cell_type_colors <- setNames(cell_colors, levels(annotation_col_df[[cell_type_title]]))
    annotation_list[[cell_type_title]] <- annotation_col_df[[cell_type_title]]
    color_list[[cell_type_title]] <- cell_type_colors
  }
  
  # // Create column annotation only if there are annotations to show
  ha_col <- NULL
  if (length(annotation_list) > 0) {
    ha_col <- do.call(HeatmapAnnotation, c(annotation_list, list(col = color_list)))
  }
  
  # // Set column split factor with custom title
  column_split_factor <- NULL
  if (show_time_annotation && time_point_title %in% colnames(annotation_col_df)) {
    column_split_factor <- annotation_col_df[[time_point_title]]
  }
  
  return(list(
    exp_mat_ordered = exp_mat_ordered,
    percent_mat_ordered = percent_mat_ordered,
    col_annotation = ha_col,
    col_split_factor = column_split_factor,
    annotation_df = annotation_col_df
  ))
}

#' Get colors with simplified approach
#' @param n_colors Number of colors needed
#' @param palette_input Either palette name (string) or color vector
#' @return Vector of colors
#' @keywords internal
get_colors_simplified <- function(n_colors, palette_input) {
  if (is.character(palette_input) && length(palette_input) == 1) {
    # Traditional palette name - use RColorBrewer
    max_colors <- min(n_colors, RColorBrewer::brewer.pal.info[palette_input, "maxcolors"])
    if (n_colors <= max_colors) {
      colors <- RColorBrewer::brewer.pal(max(3, n_colors), palette_input)[1:n_colors]
    } else {
      base_colors <- RColorBrewer::brewer.pal(max_colors, palette_input)
      colors <- rep(base_colors, length.out = n_colors)
    }
  } else if (is.character(palette_input) && length(palette_input) > 1) {
    # User provided color vector directly (from any source)
    colors <- palette_input[1:min(n_colors, length(palette_input))]
    if (length(colors) < n_colors) {
      colors <- rep(colors, length.out = n_colors)
    }
  } else {
    # Fallback to default
    colors <- RColorBrewer::brewer.pal(min(n_colors, 9), "Set1")
  }
  
  return(colors[1:n_colors])
}
