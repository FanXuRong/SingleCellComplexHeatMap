# Suppress R CMD check notes for dplyr column names
utils::globalVariables(c("GeneGroup", "Gene"))

#' Create Gene Group Annotations for Heatmap Rows
#'
#' Creates gene grouping annotations and reorders expression matrices based on gene classifications.
#'
#' @param exp_mat Expression matrix with genes as rows
#' @param percent_mat Percentage matrix with genes as rows  
#' @param gene_classification Named list where names are group labels and values are character vectors of gene names
#' @param color_palette Character string specifying palette name OR character vector of colors (default: "Set1")
#' @param sort_within_groups Logical indicating whether to sort genes within each group (default: TRUE)
#' @param annotation_title Character string for annotation title (default: "Gene Group")
#'
#' @return A list containing exp_mat_ordered (reordered expression matrix), percent_mat_ordered (reordered percentage matrix), row_annotation (ComplexHeatmap row annotation object), row_split_factor (factor for row splitting), and annotation_df (data frame with gene annotations).
#'
#' @export
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils globalVariables
#'
#' @examples
#' # Load a small example Seurat object
#' data("pbmc_small", package = "SeuratObject")
#' features <- c("CD3D", "CD79A", "MS4A1", "GZMK", "CCL5")
#' 
#' # Prepare expression matrices first
#' matrices <- prepare_expression_matrices(pbmc_small, features, group_by = "RNA_snn_res.0.8")
#' 
#' # Define gene groups
#' gene_groups <- list(
#'   "T-cell Markers" = c("CD3D", "GZMK", "CCL5"),
#'   "B-cell Markers" = c("CD79A", "MS4A1")
#' )
#' 
#' # Create gene annotations
#' annotations <- create_gene_annotations(
#'   exp_mat = matrices$exp_mat,
#'   percent_mat = matrices$percent_mat,
#'   gene_classification = gene_groups,
#'   color_palette = "Set1"
#' )
#' 
#' # Access results
#' ordered_exp_mat <- annotations$exp_mat_ordered
#' 
#' @seealso \code{\link{create_single_cell_complex_heatmap}}, \code{\link{prepare_expression_matrices}}
create_gene_annotations <- function(exp_mat, 
                                  percent_mat, 
                                  gene_classification,
                                  color_palette = "Set1",
                                  sort_within_groups = TRUE,
                                  annotation_title = "Gene Group") {
  
  # Create gene group annotations
  gene_groups_list <- lapply(names(gene_classification), function(group_name) {
    genes_in_group <- gene_classification[[group_name]]
    valid_genes <- genes_in_group[genes_in_group %in% rownames(exp_mat)]
    
    if (length(valid_genes) > 0) {
      data.frame(
        Gene = valid_genes, 
        GeneGroup = factor(rep(group_name, length(valid_genes)), 
                          levels = names(gene_classification))
      )
    } else {
      NULL
    }
  })
  
  annotation_row_df <- do.call(rbind, gene_groups_list)
  annotation_row_df$Gene <- as.character(annotation_row_df$Gene)
  
  # Order genes based on groups
  ordered_genes <- rownames(exp_mat)[rownames(exp_mat) %in% annotation_row_df$Gene]
  annotation_row_df <- annotation_row_df[match(ordered_genes, annotation_row_df$Gene), ]
  
  if (sort_within_groups) {
    annotation_row_df <- annotation_row_df %>% arrange(GeneGroup, Gene)
  } else {
    annotation_row_df <- annotation_row_df %>% arrange(GeneGroup)
  }
  
  # Reorder matrices
  exp_mat_ordered <- exp_mat[annotation_row_df$Gene, , drop = FALSE]
  percent_mat_ordered <- percent_mat[annotation_row_df$Gene, , drop = FALSE]
  
  # Enhanced color scheme handling - simplified approach
  n_groups <- length(levels(annotation_row_df$GeneGroup))
  
  if (is.character(color_palette) && length(color_palette) == 1) {
    # Traditional palette name - use RColorBrewer
    group_colors <- RColorBrewer::brewer.pal(min(n_groups, 9), color_palette)
  } else if (is.character(color_palette) && length(color_palette) > 1) {
    # User provided color vector directly (from any source: RColorBrewer, ggsci, viridis, etc.)
    group_colors <- color_palette[1:min(n_groups, length(color_palette))]
    if (length(group_colors) < n_groups) {
      group_colors <- rep(group_colors, length.out = n_groups)
    }
  } else {
    # Fallback
    group_colors <- RColorBrewer::brewer.pal(min(n_groups, 9), "Set1")
  }
  
  group_colors <- setNames(group_colors[1:n_groups], levels(annotation_row_df$GeneGroup))
  
  # Create row annotation with custom title
  annotation_args <- list()
  annotation_args[[annotation_title]] <- annotation_row_df$GeneGroup
  color_args <- list()
  color_args[[annotation_title]] <- group_colors
  
  ha_row <- do.call(rowAnnotation, c(annotation_args, list(col = color_args)))
  
  row_split_factor <- annotation_row_df$GeneGroup
  
  return(list(
    exp_mat_ordered = exp_mat_ordered,
    percent_mat_ordered = percent_mat_ordered,
    row_annotation = ha_row,
    row_split_factor = row_split_factor,
    annotation_df = annotation_row_df
  ))
}
