# Suppress R CMD check notes for dplyr/tidyr column names
utils::globalVariables(c("pct.exp", "avg.exp", "id", "avg.exp.scaled", "features.plot"))

#' Prepare Expression and Percentage Matrices from Seurat DotPlot
#'
#' Extracts and reshapes expression data from a Seurat DotPlot object into
#' matrices suitable for complex heatmap visualization.
#'
#' @param seurat_object A Seurat object containing single cell data
#' @param features Character vector of gene names to plot
#' @param group_by Character string specifying the metadata column to group by (default: "seurat_clusters")
#' @param idents Numeric or character vector specifying which cell groups to include (default: NULL for all)
#' @param split_pattern Character string used to split column names for parsing (default: "_")
#' @param time_position Integer indicating position of time point in split names (default: 1)
#' @param celltype_start Integer indicating starting position of cell type in split names (default: 2)
#'
#' @return A list containing exp_mat (matrix of scaled expression values), percent_mat (matrix of expression percentages), and dotplot_data (original DotPlot data frame).
#'
#' @export
#' @importFrom Seurat DotPlot
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom utils globalVariables
#'
#' @examples
#' # Load a small example Seurat object
#' data("pbmc_small", package = "SeuratObject")
#' features <- c("CD3D", "CD79A", "MS4A1")
#' 
#' # Basic usage
#' matrices <- prepare_expression_matrices(
#'   seurat_object = pbmc_small,
#'   features = features,
#'   group_by = "RNA_snn_res.0.8"
#' )
#' 
#' # Access the results
#' expression_matrix <- matrices$exp_mat
#' percentage_matrix <- matrices$percent_mat
#' 
#' @seealso \code{\link{create_single_cell_complex_heatmap}}
prepare_expression_matrices <- function(seurat_object, 
                                      features, 
                                      group_by = "seurat_clusters",
                                      idents = NULL,
                                      split_pattern = "_",
                                      time_position = 1,
                                      celltype_start = 2) {
  
  # Parameter validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }
  
  if (!group_by %in% colnames(seurat_object@meta.data)) {
    stop("group_by column '", group_by, "' not found in metadata")
  }
  
  # Add error handling for DotPlot
  tryCatch({
    if (is.null(idents)) {
      p <- DotPlot(object = seurat_object, features = features, group.by = group_by)
    } else {
      p <- DotPlot(object = seurat_object, features = features, group.by = group_by, idents = idents)
    }
  }, error = function(e) {
    stop("Error creating DotPlot: ", e$message)
  })
  
  df <- p$data
  
  # Create expression matrix
  exp_mat <- df %>%
    dplyr::select(-pct.exp, -avg.exp) %>%
    tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  
  rownames(exp_mat) <- exp_mat$features.plot
  exp_mat <- exp_mat[, -1] %>% as.matrix()
  
  # Create percentage matrix
  percent_mat <- df %>%
    dplyr::select(-avg.exp, -avg.exp.scaled) %>%
    tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()
  
  rownames(percent_mat) <- percent_mat$features.plot
  percent_mat <- percent_mat[, -1] %>% as.matrix()
  
  return(list(
    exp_mat = exp_mat,
    percent_mat = percent_mat,
    dotplot_data = df
  ))
}
