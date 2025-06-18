# Suppress R CMD check notes for dplyr column names
utils::globalVariables(c("TimePoint", "CellType"))

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
#' @param time_color_palette Character string specifying RColorBrewer palette for time points (default: "Accent")
#' @param celltype_color_palette Character string specifying RColorBrewer palette for cell types (default: "Dark2")
#' @param show_time_annotation Logical indicating whether to show time point annotation (default: TRUE)
#' @param show_celltype_annotation Logical indicating whether to show cell type annotation (default: TRUE)
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
#' \dontrun{
#' # Prepare expression matrices first
#' matrices <- prepare_expression_matrices(seurat_obj, features, group_by = "timepoint_celltype")
#' 
#' # Create cell annotations with custom ordering
#' col_annotations <- create_cell_annotations(
#'   exp_mat = matrices$exp_mat,
#'   percent_mat = matrices$percent_mat,
#'   split_pattern = "_",
#'   time_points_order = c("0h", "6h", "12h", "24h"),
#'   cell_types_order = c("T_cell", "B_cell", "Monocyte"),
#'   time_color_palette = "Set2",
#'   celltype_color_palette = "Dark2"
#' )
#' 
#' # Access results
#' ordered_exp_mat <- col_annotations$exp_mat_ordered
#' col_annotation <- col_annotations$col_annotation
#' col_split_factor <- col_annotations$col_split_factor
#' }
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
                                  show_celltype_annotation = TRUE) {
  
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
          # If no underscore found, use the whole string as cell type
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
    
    # // Create annotation data frame
    annotation_col_df <- data.frame(id = col_ids)
    
    if (show_time_annotation) {
      annotation_col_df$TimePoint <- factor(time_points_extracted, levels = time_points_order)
    }
    
    if (show_celltype_annotation) {
      annotation_col_df$CellType <- factor(cell_types_extracted, levels = cell_types_order)
    }
    
    # // Arrange data frame based on what annotations are shown
    if (show_time_annotation && show_celltype_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(TimePoint, CellType)
    } else if (show_time_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(TimePoint)
    } else if (show_celltype_annotation) {
      annotation_col_df <- annotation_col_df %>% arrange(CellType)
    }
  } else {
    # // If no annotations are shown, just keep original order
    annotation_col_df <- data.frame(id = col_ids)
  }
  
  # // Reorder matrices
  exp_mat_ordered <- exp_mat[, annotation_col_df$id, drop = FALSE]
  percent_mat_ordered <- percent_mat[, annotation_col_df$id, drop = FALSE]
  
  # // Create color schemes and annotations
  annotation_list <- list()
  color_list <- list()
  
  if (show_time_annotation && "TimePoint" %in% colnames(annotation_col_df)) {
    n_timepoints <- length(levels(annotation_col_df$TimePoint))
    max_colors <- min(n_timepoints, RColorBrewer::brewer.pal.info[time_color_palette, "maxcolors"])
    if (n_timepoints <= max_colors) {
      time_colors <- RColorBrewer::brewer.pal(max(3, n_timepoints), time_color_palette)[1:n_timepoints]
    } else {
      base_colors <- RColorBrewer::brewer.pal(max_colors, time_color_palette)
      time_colors <- rep(base_colors, length.out = n_timepoints)
    }
    
    time_point_colors <- setNames(time_colors, levels(annotation_col_df$TimePoint))
    annotation_list$TimePoint <- annotation_col_df$TimePoint
    color_list$TimePoint <- time_point_colors
  }
  
  if (show_celltype_annotation && "CellType" %in% colnames(annotation_col_df)) {
    n_celltypes <- length(levels(annotation_col_df$CellType))
    max_colors <- min(n_celltypes, RColorBrewer::brewer.pal.info[celltype_color_palette, "maxcolors"])
    if (n_celltypes <= max_colors) {
      cell_colors <- RColorBrewer::brewer.pal(max(3, n_celltypes), celltype_color_palette)[1:n_celltypes]
    } else {
      base_colors <- RColorBrewer::brewer.pal(max_colors, celltype_color_palette)
      cell_colors <- rep(base_colors, length.out = n_celltypes)
    }
    
    cell_type_colors <- setNames(cell_colors, levels(annotation_col_df$CellType))
    annotation_list$CellType <- annotation_col_df$CellType
    color_list$CellType <- cell_type_colors
  }
  
  # // Create column annotation only if there are annotations to show
  ha_col <- NULL
  if (length(annotation_list) > 0) {
    ha_col <- do.call(HeatmapAnnotation, c(annotation_list, list(col = color_list)))
  }
  
  # // Set column split factor
  column_split_factor <- NULL
  if (show_time_annotation && "TimePoint" %in% colnames(annotation_col_df)) {
    column_split_factor <- annotation_col_df$TimePoint
  }
  
  return(list(
    exp_mat_ordered = exp_mat_ordered,
    percent_mat_ordered = percent_mat_ordered,
    col_annotation = ha_col,
    col_split_factor = column_split_factor,
    annotation_df = annotation_col_df
  ))
}
