% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_gene_annotations.R
\name{create_gene_annotations}
\alias{create_gene_annotations}
\title{Create Gene Group Annotations for Heatmap Rows}
\usage{
create_gene_annotations(
  exp_mat,
  percent_mat,
  gene_classification,
  color_palette = "Set1",
  sort_within_groups = TRUE,
  annotation_title = "Gene Group"
)
}
\arguments{
\item{exp_mat}{Expression matrix with genes as rows}

\item{percent_mat}{Percentage matrix with genes as rows}

\item{gene_classification}{Named list where names are group labels and values are character vectors of gene names}

\item{color_palette}{Character string specifying palette name OR character vector of colors (default: "Set1")}

\item{sort_within_groups}{Logical indicating whether to sort genes within each group (default: TRUE)}

\item{annotation_title}{Character string for annotation title (default: "Gene Group")}
}
\value{
A list containing exp_mat_ordered (reordered expression matrix), percent_mat_ordered (reordered percentage matrix), row_annotation (ComplexHeatmap row annotation object), row_split_factor (factor for row splitting), and annotation_df (data frame with gene annotations).
}
\description{
Creates gene grouping annotations and reorders expression matrices based on gene classifications.
}
\examples{
# Load a small example Seurat object
data("pbmc_small", package = "SeuratObject")
features <- c("CD3D", "CD79A", "MS4A1", "GZMK", "CCL5")

# Prepare expression matrices first
matrices <- prepare_expression_matrices(pbmc_small, features, group_by = "RNA_snn_res.0.8")

# Define gene groups
gene_groups <- list(
  "T-cell Markers" = c("CD3D", "GZMK", "CCL5"),
  "B-cell Markers" = c("CD79A", "MS4A1")
)

# Create gene annotations
annotations <- create_gene_annotations(
  exp_mat = matrices$exp_mat,
  percent_mat = matrices$percent_mat,
  gene_classification = gene_groups,
  color_palette = "Set1"
)

# Access results
ordered_exp_mat <- annotations$exp_mat_ordered

}
\seealso{
\code{\link{create_single_cell_complex_heatmap}}, \code{\link{prepare_expression_matrices}}
}
