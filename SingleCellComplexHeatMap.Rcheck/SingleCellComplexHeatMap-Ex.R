pkgname <- "SingleCellComplexHeatMap"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SingleCellComplexHeatMap-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SingleCellComplexHeatMap')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("create_cell_annotations")
### * create_cell_annotations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_cell_annotations
### Title: Create Cell Type and Time Point Annotations for Heatmap Columns
### Aliases: create_cell_annotations

### ** Examples

## Not run: 
##D # Prepare expression matrices first
##D matrices <- prepare_expression_matrices(seurat_obj, features, group_by = "timepoint_celltype")
##D 
##D # Create cell annotations with custom ordering
##D col_annotations <- create_cell_annotations(
##D   exp_mat = matrices$exp_mat,
##D   percent_mat = matrices$percent_mat,
##D   split_pattern = "_",
##D   time_points_order = c("0h", "6h", "12h", "24h"),
##D   cell_types_order = c("T_cell", "B_cell", "Monocyte"),
##D   time_color_palette = "Set2",
##D   celltype_color_palette = "Dark2"
##D )
##D 
##D # Access results
##D ordered_exp_mat <- col_annotations$exp_mat_ordered
##D col_annotation <- col_annotations$col_annotation
##D col_split_factor <- col_annotations$col_split_factor
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_cell_annotations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_gene_annotations")
### * create_gene_annotations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_gene_annotations
### Title: Create Gene Group Annotations for Heatmap Rows
### Aliases: create_gene_annotations

### ** Examples

## Not run: 
##D # Prepare expression matrices first
##D matrices <- prepare_expression_matrices(seurat_obj, features, group_by = "cell_type")
##D 
##D # Define gene groups
##D gene_groups <- list(
##D   "Immune_Response" = c("Gene1", "Gene2"),
##D   "Cell_Cycle" = c("Gene3", "Gene4"),
##D   "Metabolism" = c("Gene5", "Gene6")
##D )
##D 
##D # Create gene annotations
##D annotations <- create_gene_annotations(
##D   exp_mat = matrices$exp_mat,
##D   percent_mat = matrices$percent_mat,
##D   gene_classification = gene_groups,
##D   color_palette = "Set1"
##D )
##D 
##D # Access results
##D ordered_exp_mat <- annotations$exp_mat_ordered
##D ordered_percent_mat <- annotations$percent_mat_ordered
##D row_annotation <- annotations$row_annotation
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_gene_annotations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prepare_expression_matrices")
### * prepare_expression_matrices

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prepare_expression_matrices
### Title: Prepare Expression and Percentage Matrices from Seurat DotPlot
### Aliases: prepare_expression_matrices

### ** Examples

## Not run: 
##D # Basic usage
##D matrices <- prepare_expression_matrices(
##D   seurat_object = my_seurat,
##D   features = c("Gene1", "Gene2", "Gene3"),
##D   group_by = "cell_type"
##D )
##D 
##D # Advanced usage with specific cell groups
##D matrices <- prepare_expression_matrices(
##D   seurat_object = my_seurat,
##D   features = gene_list,
##D   group_by = "timepoint_celltype",
##D   idents = c("0h_T_cell", "6h_T_cell", "12h_T_cell"),
##D   split_pattern = "_"
##D )
##D 
##D # Access the results
##D expression_matrix <- matrices$exp_mat
##D percentage_matrix <- matrices$percent_mat
##D original_data <- matrices$dotplot_data
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prepare_expression_matrices", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
