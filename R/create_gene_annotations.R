# Suppress R CMD check notes for dplyr column names
utils::globalVariables(c("GeneGroup", "Gene"))

#' Create Gene Group Annotations for Heatmap Rows
#'
#' Creates gene grouping annotations and reorders expression matrices based on gene classifications.
#'
#' @param exp_mat Expression matrix with genes as rows
#' @param percent_mat Percentage matrix with genes as rows  
#' @param gene_classification Named list where names are group labels and values are character vectors of gene names
#' @param color_palette Character string specifying RColorBrewer palette name (default: "Set1")
#' @param sort_within_groups Logical indicating whether to sort genes within each group (default: TRUE)
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
#' \dontrun{
#' # Prepare expression matrices first
#' matrices <- prepare_expression_matrices(seurat_obj, features, group_by = "cell_type")
#' 
#' # Define gene groups
#' gene_groups <- list(
#'   "Immune_Response" = c("Gene1", "Gene2"),
#'   "Cell_Cycle" = c("Gene3", "Gene4"),
#'   "Metabolism" = c("Gene5", "Gene6")
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
#' ordered_percent_mat <- annotations$percent_mat_ordered
#' row_annotation <- annotations$row_annotation
#' }
#' 
#' @seealso \code{\link{create_single_cell_complex_heatmap}}, \code{\link{prepare_expression_matrices}}
create_gene_annotations <- function(exp_mat, 
                                  percent_mat, 
                                  gene_classification,
                                  color_palette = "Set1",
                                  sort_within_groups = TRUE) {
  
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
  
  # Create color scheme
  n_groups <- length(levels(annotation_row_df$GeneGroup))
  group_colors <- setNames(
    RColorBrewer::brewer.pal(min(n_groups, 9), color_palette),
    levels(annotation_row_df$GeneGroup)
  )
  
  # Create row annotation
  ha_row <- rowAnnotation(
    GeneGroup = annotation_row_df$GeneGroup,
    col = list(GeneGroup = group_colors)
  )
  
  row_split_factor <- annotation_row_df$GeneGroup
  
  return(list(
    exp_mat_ordered = exp_mat_ordered,
    percent_mat_ordered = percent_mat_ordered,
    row_annotation = ha_row,
    row_split_factor = row_split_factor,
    annotation_df = annotation_row_df
  ))
}
