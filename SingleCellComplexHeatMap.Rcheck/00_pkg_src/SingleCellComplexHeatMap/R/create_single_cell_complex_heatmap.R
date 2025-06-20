#' Create Complex Heatmap for Single Cell Expression Data
#'
#' @description
#' Creates a complex heatmap that displays both gene expression levels (as color intensity) 
#' and expression percentages (as circle sizes) for single cell RNA-seq data.
#' This function provides extensive customization options while maintaining ease of use.
#'
#' @param seurat_object A Seurat object containing single cell data
#' @param features Character vector of gene names to plot
#' @param gene_classification Named list where names are group labels and values are character vectors of gene names (default: NULL for no gene grouping)
#' @param group_by Character string specifying the metadata column to group by (default: "seurat_clusters")
#' @param idents Numeric or character vector specifying which cell groups to include (default: NULL for all)
#' @param time_points_order Character vector specifying order of time points. Only affects display order, not data filtering (default: NULL for automatic)
#' @param cell_types_order Character vector specifying order of cell types. Only affects display order, not data filtering (default: NULL for automatic)
#' @param color_range Numeric vector specifying color mapping break points for expression values. Its length must match color_palette if color_palette is a vector. (default: c(-1, 0, 2))
#' @param color_palette Character vector specifying colors for expression heatmap. Its length must match color_range. If NULL, a default palette (viridis or color_palette_main) is generated to match color_range length (default: NULL)
#' @param max_circle_size Numeric specifying maximum circle radius in mm. This applies to the highest percentage value in percentage_breaks (default: 2)
#' @param row_fontsize Numeric specifying row name font size (default: 8)
#' @param col_fontsize Numeric specifying column name font size (default: 9)
#' @param col_name_rotation Numeric specifying column name rotation angle (default: 90)
#' @param row_title_fontsize Numeric specifying row title font size (default: 10)
#' @param col_title_fontsize Numeric specifying column title font size (default: 10)
#' @param show_heatmap_legend Logical indicating whether to show heatmap legend (default: TRUE)
#' @param show_percentage_legend Logical indicating whether to show percentage legend (default: TRUE)
#' @param legend_side Character string specifying legend position (default: "right")
#' @param cell_border_color Character string specifying cell border color (default: "grey80")
#' @param split_pattern Character string used to split column names for parsing (default: "_")
#' @param gene_color_palette Character string specifying RColorBrewer palette for gene groups (default: "Set1")
#' @param time_color_palette Character string specifying RColorBrewer palette for time points (default: "Accent")
#' @param celltype_color_palette Character string specifying RColorBrewer palette for cell types (default: "Dark2")
#' @param show_gene_grouping Logical indicating whether to show gene grouping (default: TRUE if gene_classification provided)
#' @param show_time_annotation Logical indicating whether to show time point annotation (default: TRUE)
#' @param show_celltype_annotation Logical indicating whether to show cell type annotation (default: TRUE)
#' @param split_by Character string specifying how to split columns: "time", "celltype", or "none" (default: "time")
#' @param merge_legends Logical indicating whether to merge legends (default: TRUE)
#' @param percentage_legend_title Character string for percentage legend title (default: "Expression %")
#' @param percentage_legend_labels Character vector for percentage legend labels
#' @param percentage_breaks Numeric vector specifying actual percentage values corresponding to labels
#' @param return_data Logical; if TRUE, return underlying data instead of drawing only
#' @param save_plot File path to save the drawn heatmap (PNG)
#' @param plot_width Numeric; width in inches for saving
#' @param plot_height Numeric; height in inches for saving
#' @param plot_dpi Numeric; resolution (DPI) for saved plot
#' @param assay Seurat assay name to extract data from
#' @param slot Seurat slot name within assay (e.g., "scale.data", "data")
#' @param cluster_cells Logical; whether to cluster columns (cells)
#' @param cluster_features Logical; whether to cluster rows (features)
#' @param clustering_distance_rows Distance metric for row clustering
#' @param clustering_distance_cols Distance metric for column clustering
#' @param clustering_method_rows Clustering method for rows
#' @param clustering_method_cols Clustering method for columns
#' @param color_palette_main Fallback color palette when viridis unavailable
#' @param annotation_colors Named list of custom annotation colors
#' @param show_feature_names Logical; whether to show feature (row) names
#' @param feature_names_gp gpar object controlling feature name appearance
#' @param legend_title Character; title for main heatmap legend
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap()
#'
#' @return A ComplexHeatmap object. If return_data is TRUE, returns a list containing the heatmap object and underlying data matrices.
#' @export
#' @importFrom ComplexHeatmap Heatmap Legend draw pindex
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.rect grid.circle gpar unit
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom grDevices png dev.off
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils globalVariables
create_single_cell_complex_heatmap <- function(seurat_object,
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
                                             percentage_legend_title = "Expression %",
                                             percentage_legend_labels = c("0%", "25%", "50%", "75%", "100%"),
                                             percentage_breaks = NULL,
                                             return_data = FALSE,
                                             save_plot = NULL,
                                             plot_width = 10,
                                             plot_height = 8,
                                             plot_dpi = 300,
                                             # Data Source Control
                                             assay = NULL,
                                             slot = 'scale.data',
                                             # Clustering Control
                                             cluster_cells = TRUE,
                                             cluster_features = TRUE,
                                             clustering_distance_rows = "euclidean",
                                             clustering_distance_cols = "euclidean", 
                                             clustering_method_rows = "complete",
                                             clustering_method_cols = "complete",
                                             # Color Mapping Control
                                             color_palette_main = c("blue", "white", "red"),
                                             annotation_colors = NULL,
                                             # Label and Legend Control
                                             show_feature_names = TRUE,
                                             feature_names_gp = NULL,
                                             legend_title = "Expression",
                                             # Ultimate flexibility
                                             ...) {
  
  # Parameter validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }
  
  if (!is.character(features) || length(features) == 0) {
    stop("features must be a non-empty character vector")
  }
  
  # Check if features exist in the Seurat object
  missing_features <- features[!features %in% rownames(seurat_object)]
  if (length(missing_features) > 0) {
    warning("The following features are not found in the Seurat object: ", 
            paste(missing_features, collapse = ", "))
    features <- features[features %in% rownames(seurat_object)]
    if (length(features) == 0) {
      stop("No valid features found in the Seurat object")
    }
  }
  
  # Validate color_range
  if (!is.numeric(color_range) || length(color_range) < 3) {
    stop("color_range must be a numeric vector with at least 3 values")
  }
  
  # Validate group_by column exists
  if (!group_by %in% colnames(seurat_object@meta.data)) {
    stop("group_by column '", group_by, "' not found in seurat_object metadata")
  }
  
  # Set default assay if not provided
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_object)
  }
  
  # Validate assay exists
  if (!assay %in% names(seurat_object@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }
  
  # Set default for show_gene_grouping based on gene_classification
  if (is.null(show_gene_grouping)) {
    show_gene_grouping <- !is.null(gene_classification)
  }
  
  # Set default color palette if not provided
  if (is.null(color_palette)) {
    if (requireNamespace("viridis", quietly = TRUE)) {
      color_palette <- viridis::viridis(20)[c(1, 10, 20)]
    } else {
      color_palette <- color_palette_main
    }
  }
  
  # Set default feature names gpar
  if (is.null(feature_names_gp)) {
    feature_names_gp <- gpar(fontsize = row_fontsize)
  }
  
  # Determine the actual color palette for the expression heatmap
  actual_expression_palette <- color_palette
  if (is.null(actual_expression_palette)) {
    if (requireNamespace("viridis", quietly = TRUE)) {
      # Generate a viridis palette with the same number of colors as breaks in color_range
      actual_expression_palette <- viridis::viridis(length(color_range))
    } else {
      # Fallback to color_palette_main
      actual_expression_palette <- color_palette_main
    }
  }

  # Validate that color_range and actual_expression_palette are compatible in length
  if (length(color_range) != length(actual_expression_palette)) {
    stop(paste0("Length of color_range (", length(color_range),
                ") must match the length of the determined color palette for expression (", length(actual_expression_palette), ").\n",
                "Current color_range: [", paste(color_range, collapse=", "), "].\n",
                "Current color_palette for expression has ", length(actual_expression_palette), " colors: [", paste(actual_expression_palette, collapse=", "), "].\n",
                "Please ensure color_range and color_palette (if provided) have the same number of elements, ",
                "or if using defaults, ensure color_palette_main's length matches color_range if viridis is unavailable."))
  }

  # Parse percentage breaks from labels if not provided
  if (is.null(percentage_breaks)) {
    # Extract numeric values from percentage labels
    numeric_labels <- gsub("[^0-9.]", "", percentage_legend_labels)
    percentage_breaks <- as.numeric(numeric_labels)
    if(any(is.na(percentage_breaks))) {
      stop("Could not parse numeric values from percentage_legend_labels. Ensure labels are like '25%' or '25', or set percentage_breaks manually.")
    }
  }
  
  # Validate that breaks and labels have same length
  if (length(percentage_breaks) != length(percentage_legend_labels)) {
    stop("percentage_breaks and percentage_legend_labels must have the same length")
  }
  
  # Ensure breaks are in ascending order
  if (!all(diff(percentage_breaks) >= 0)) {
    warning("percentage_breaks should be in ascending order for proper legend display")
  }
  
  # Validate that color_range and color_palette are compatible
  if (length(color_range) != length(color_palette)) {
    # Auto-adjust color_palette to match color_range if needed
    if (length(color_palette) == 3 && length(color_range) > 3) {
      # Extend palette for longer color_range
      warning("color_palette has fewer colors than color_range. Interpolating additional colors.")
    } else if (length(color_palette) > length(color_range)) {
      # Trim palette to match color_range
      color_palette <- color_palette[1:length(color_range)]
      warning("color_palette trimmed to match color_range length")
    }
  }
  
  # Step 1: Prepare expression matrices using DotPlot approach for compatibility
  matrices <- prepare_expression_matrices(
    seurat_object = seurat_object,
    features = features,
    group_by = group_by,
    idents = idents,
    split_pattern = split_pattern
  )
  
  # Alternative: Extract data directly from assay if user specifies different slot
  if (slot != 'scale.data') {
    tryCatch({
      direct_data <- GetAssayData(seurat_object, assay = assay, slot = slot)
      if (all(features %in% rownames(direct_data))) {
        # Extract relevant cells based on group_by and idents
        if (!is.null(idents)) {
          relevant_cells <- rownames(seurat_object@meta.data)[seurat_object@meta.data[[group_by]] %in% idents]
        } else {
          relevant_cells <- colnames(seurat_object)
        }
        
        # Use the direct data for expression matrix
        exp_mat_current <- as.matrix(direct_data[features, relevant_cells, drop = FALSE])
        # Keep percentage matrix from DotPlot for circle sizes
        percent_mat_current <- matrices$percent_mat
      } else {
        exp_mat_current <- matrices$exp_mat
        percent_mat_current <- matrices$percent_mat
        warning("Some features not found in specified slot, using DotPlot data")
      }
    }, error = function(e) {
      warning("Could not extract data from specified slot, using DotPlot data: ", e$message)
      exp_mat_current <- matrices$exp_mat
      percent_mat_current <- matrices$percent_mat
    })
  } else {
    exp_mat_current <- matrices$exp_mat
    percent_mat_current <- matrices$percent_mat
  }
  
  row_annotation <- NULL
  row_split_factor <- NULL
  
  # Step 2: Create gene annotations if requested
  if (show_gene_grouping && !is.null(gene_classification)) {
    gene_annotations <- create_gene_annotations(
      exp_mat = exp_mat_current,
      percent_mat = percent_mat_current,
      gene_classification = gene_classification,
      color_palette = gene_color_palette
    )
    exp_mat_current <- gene_annotations$exp_mat_ordered
    percent_mat_current <- gene_annotations$percent_mat_ordered
    row_annotation <- gene_annotations$row_annotation
    row_split_factor <- gene_annotations$row_split_factor
    
    # Apply custom annotation colors if provided
    if (!is.null(annotation_colors) && "GeneGroup" %in% names(annotation_colors)) {
      # Update row annotation with custom colors
      gene_annotations$row_annotation@anno_list$GeneGroup@color_mapping@colors <- annotation_colors$GeneGroup
    }
  }
  
  col_annotation <- NULL
  col_split_factor <- NULL
  
  # Step 3: Create cell annotations if requested
  if (show_time_annotation || show_celltype_annotation) {
    cell_annotations <- create_cell_annotations(
      exp_mat = exp_mat_current,
      percent_mat = percent_mat_current,
      split_pattern = split_pattern,
      time_points_order = time_points_order,
      cell_types_order = cell_types_order,
      time_color_palette = time_color_palette,
      celltype_color_palette = celltype_color_palette,
      show_time_annotation = show_time_annotation,
      show_celltype_annotation = show_celltype_annotation
    )
    exp_mat_current <- cell_annotations$exp_mat_ordered
    percent_mat_current <- cell_annotations$percent_mat_ordered
    col_annotation <- cell_annotations$col_annotation
    
    # Apply custom annotation colors if provided
    if (!is.null(annotation_colors)) {
      if ("TimePoint" %in% names(annotation_colors) && show_time_annotation) {
        # Update time point colors
        if (!is.null(col_annotation)) {
          col_annotation@anno_list$TimePoint@color_mapping@colors <- annotation_colors$TimePoint
        }
      }
      if ("CellType" %in% names(annotation_colors) && show_celltype_annotation) {
        # Update cell type colors
        if (!is.null(col_annotation)) {
          col_annotation@anno_list$CellType@color_mapping@colors <- annotation_colors$CellType
        }
      }
    }
    
    # Set column split based on split_by parameter
    if (split_by == "time" && show_time_annotation) {
      col_split_factor <- cell_annotations$col_split_factor
    } else if (split_by == "celltype" && show_celltype_annotation) {
      col_split_factor <- cell_annotations$annotation_df$CellType
    } else {
      col_split_factor <- NULL
    }
  }
  
  # Step 4: Create color function for expression heatmap
  col_fun <- circlize::colorRamp2(color_range, actual_expression_palette)
  
  # Step 5: Define cell and layer functions with custom percentage mapping
  cell_fun <- function(j, i, x, y, w, h, fill) {
    if (!is.na(cell_border_color)) {
      grid.rect(x = x, y = y, width = w, height = h,
               gp = gpar(col = cell_border_color, fill = NA))
    }
  }
  
  layer_fun <- function(j, i, x, y, w, h, fill) {
    percentage_values_slice <- pindex(percent_mat_current, i, j)
    expression_values_slice <- pindex(exp_mat_current, i, j)
    
    valid_indices <- !is.na(percentage_values_slice) & !is.na(expression_values_slice)
    
    if (any(valid_indices)) {
      # 强制转换为数值型，防止字符型导致报错
      percentage_numeric <- as.numeric(percentage_values_slice[valid_indices])
      # Map percentage values to circle sizes based on custom breaks
      max_break <- max(percentage_breaks)
      min_break <- min(percentage_breaks)
      
      # Normalize percentage values to 0-1 scale based on breaks
      normalized_percentages <- pmax(0, pmin(1, 
        (percentage_numeric - min_break) / (max_break - min_break)
      ))
      
      grid.circle(
        x = x[valid_indices],
        y = y[valid_indices],
        r = normalized_percentages * unit(max_circle_size, "mm"),
        gp = gpar(fill = col_fun(expression_values_slice[valid_indices]), col = NA)
      )
    }
  }
  
  # Step 6: Create percentage legend if requested
  legend_list <- list()
  if (show_percentage_legend) {
    # Use custom breaks for legend
    max_break <- max(percentage_breaks)
    min_break <- min(percentage_breaks)
    
    # Calculate normalized values for legend circles
    legend_values <- (percentage_breaks - min_break) / (max_break - min_break)
    
    legend_graphics <- lapply(legend_values, function(val) {
      function(x, y, w, h) grid.circle(x = x, y = y, 
                                      r = val * unit(max_circle_size, "mm"),
                                      gp = gpar(fill = "black"))
    })
    
    legend_list <- list(
      Legend(
        labels = percentage_legend_labels,
        title = percentage_legend_title,
        graphics = legend_graphics
      )
    )
  }
  
  # Step 7: Prepare clustering parameters
  cluster_rows_param <- if (show_gene_grouping && !is.null(gene_classification)) {
    FALSE  # Don't cluster if genes are grouped
  } else {
    cluster_features
  }
  
  cluster_columns_param <- cluster_cells
  
  # Prepare clustering distance and method parameters
  clustering_params <- list()
  if (cluster_rows_param) {
    clustering_params$clustering_distance_rows <- clustering_distance_rows
    clustering_params$clustering_method_rows <- clustering_method_rows
  }
  if (cluster_columns_param) {
    clustering_params$clustering_distance_columns <- clustering_distance_cols
    clustering_params$clustering_method_columns <- clustering_method_cols
  }
  
  # Step 8: Create heatmap with all parameters
  base_params <- list(
    matrix = exp_mat_current,
    name = legend_title,
    col = col_fun,
    rect_gp = gpar(type = "none"),
    cell_fun = cell_fun,
    cluster_rows = cluster_rows_param,
    cluster_columns = cluster_columns_param,
    show_row_names = show_feature_names,
    row_names_gp = feature_names_gp,
    column_names_gp = gpar(fontsize = col_fontsize),
    column_names_rot = col_name_rotation,
    top_annotation = col_annotation,
    left_annotation = row_annotation,
    row_split = row_split_factor,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = row_title_fontsize, fontface = "bold"),
    column_split = col_split_factor,
    column_title_gp = gpar(fontsize = col_title_fontsize, fontface = "bold"),
    show_heatmap_legend = show_heatmap_legend,
    heatmap_legend_param = list(
      title = legend_title,
      col_fun = col_fun, # Uses the validated col_fun
      at = color_range,  # Break points match col_fun's breaks
      labels = as.character(color_range) # Labels for these breaks
    ),
    layer_fun = layer_fun
  )
  
  # Merge clustering parameters
  all_params <- c(base_params, clustering_params)
  
  # Add any additional parameters from ...
  dots <- list(...)
  all_params <- c(all_params, dots)
  
  # Create heatmap
  ht_plot <- do.call(Heatmap, all_params)
  
  # Step 9: Draw heatmap with legends
  draw(ht_plot, 
       heatmap_legend_side = legend_side, 
       annotation_legend_side = legend_side,
       merge_legends = merge_legends, 
       annotation_legend_list = legend_list)
  
  # Option to return underlying data
  if (return_data) {
    return(list(
      heatmap = ht_plot,
      expression_matrix = exp_mat_current,
      percentage_matrix = percent_mat_current,
      gene_annotations = if(exists("gene_annotations")) gene_annotations else NULL,
      cell_annotations = if(exists("cell_annotations")) cell_annotations else NULL
    ))
  }
  
  # Option to save plot
  if (!is.null(save_plot)) {
    grDevices::png(save_plot, width = plot_width, height = plot_height, 
        units = "in", res = plot_dpi)
    draw(ht_plot, 
         heatmap_legend_side = legend_side, 
         annotation_legend_side = legend_side,
         merge_legends = merge_legends, 
         annotation_legend_list = legend_list)
    grDevices::dev.off()
    message("Plot saved to: ", save_plot)
  }
  
  return(ht_plot)
}
