# SingleCellComplexHeatMap

A R package for creating complex heatmaps for single cell RNA-seq data that simultaneously display gene expression levels (as color intensity) and expression percentages (as circle sizes).

## Features

- **Dual Information Display**: Shows both expression intensity and percentage in a single visualization
- **Gene Grouping**: Organize genes by functional categories with color-coded annotations
- **Time Course Support**: Handle time series data with proper ordering and visualization
- **Cell Type Annotations**: Annotate and group different cell types
- **Flexible Display Options**: Show/hide different annotation layers (gene groups, time points, cell types)
- **Highly Customizable**: Extensive parameters for colors, sizes, fonts, and layouts
- **Seurat Integration**: Direct integration with Seurat objects
- **Enhanced Color Support**: Compatible with multiple R color packages (RColorBrewer, ggsci, viridis, custom colors)
- **Customizable Annotations**: User-defined titles for gene groups, time points, and cell types
- **Flexible Visual Control**: Optional cell borders and column annotations
- **Gene Name Mapping**: Custom display names for genes on y-axis

## Installation

```r
# Install from GitHub (when available)
# devtools::install_github("xuecheng328/SingleCellComplexHeatMap")

# For now, install locally
devtools::install_local("/path/to/SingleCellComplexHeatMap")
```

## Function Parameters

### Main Function: `create_single_cell_complex_heatmap()`

#### Required Parameters
- **`seurat_object`** - A Seurat object containing single cell data
- **`features`** - Character vector of gene names to plot

#### Data Organization Parameters
- **`gene_classification`** - Named list where names are group labels and values are character vectors of gene names (default: `NULL` for no gene grouping)
- **`group_by`** - Character string specifying the metadata column to group by (default: `"seurat_clusters"`)
- **`idents`** - Numeric or character vector specifying which cell groups to include (default: `NULL` for all)

#### Order Control Parameters
- **`time_points_order`** - Character vector specifying order of time points (default: `NULL` for automatic)
- **`cell_types_order`** - Character vector specifying order of cell types (default: `NULL` for automatic)

#### Color and Style Parameters
- **`color_range`** - Numeric vector of length 3 or 4 specifying color mapping range (default: `c(-1, 0, 2)`)
- **`color_palette`** - Function or character vector specifying colors (default: viridis palette)
- **`max_circle_size`** - Numeric specifying maximum circle radius in mm (default: `2`)

#### Font Parameters
- **`row_fontsize`** - Numeric specifying row name font size (default: `8`)
- **`col_fontsize`** - Numeric specifying column name font size (default: `9`)
- **`col_name_rotation`** - Numeric specifying column name rotation angle (default: `90`)
- **`row_title_fontsize`** - Numeric specifying row title font size (default: `10`)
- **`col_title_fontsize`** - Numeric specifying column title font size (default: `10`)

#### Legend Parameters
- **`show_heatmap_legend`** - Logical indicating whether to show heatmap legend (default: `TRUE`)
- **`show_percentage_legend`** - Logical indicating whether to show percentage legend (default: `TRUE`)
- **`legend_side`** - Character string specifying legend position (default: `"right"`)
- **`merge_legends`** - Logical indicating whether to merge legends (default: `TRUE`)
- **`percentage_legend_title`** - Character string for percentage legend title (default: `"Expression %"`)
- **`percentage_legend_labels`** - Character vector for percentage legend labels (default: `c("0%", "25%", "50%", "75%", "100%")`)

#### Visual Appearance Parameters
- **`cell_border_color`** - Character string specifying cell border color (default: `"grey80"`)
- **`show_cell_borders`** - **NEW**: Logical indicating whether to show cell border lines (default: `TRUE`)

#### Data Parsing Parameters
- **`split_pattern`** - Character string used to split column names for parsing (default: `"_"`)

#### Color Palette Parameters
- **`gene_color_palette`** - **ENHANCED**: Character string specifying palette name OR character vector of colors for gene groups (default: `"Set1"`)
- **`time_color_palette`** - **ENHANCED**: Character string specifying palette name OR character vector of colors for time points (default: `"Accent"`)
- **`celltype_color_palette`** - **ENHANCED**: Character string specifying palette name OR character vector of colors for cell types (default: `"Dark2"`)

#### Display Control Parameters
- **`show_gene_grouping`** - Logical indicating whether to show gene grouping (default: `TRUE` if gene_classification provided)
- **`show_time_annotation`** - Logical indicating whether to show time point annotation (default: `TRUE`)
- **`show_celltype_annotation`** - Logical indicating whether to show cell type annotation (default: `TRUE`)
- **`show_column_annotation`** - **NEW**: Logical indicating whether to show column annotations (default: `TRUE`)
- **`split_by`** - Character string specifying how to split columns: `"time"`, `"celltype"`, or `"none"` (default: `"time"`)

#### Customization Parameters (NEW)
- **`gene_group_title`** - Character string for gene group annotation title (default: `"Gene Group"`)
- **`time_point_title`** - Character string for time point annotation title (default: `"Time Point"`)
- **`cell_type_title`** - Character string for cell type annotation title (default: `"Cell Type"`)
- **`gene_name_mapping`** - Named character vector for mapping gene names, where names are original gene names and values are display names (default: `NULL`)

### Helper Functions

#### `prepare_expression_matrices()`
```r
prepare_expression_matrices(
  seurat_object,           # Required: Seurat object
  features,                # Required: Gene names
  group_by = "seurat_clusters",
  idents = NULL,
  split_pattern = "_",
  time_position = 1,
  celltype_start = 2
)
```

#### `create_gene_annotations()`
```r
create_gene_annotations(
  exp_mat,                 # Required: Expression matrix
  percent_mat,             # Required: Percentage matrix
  gene_classification,     # Required: Gene groups
  color_palette = "Set1",  # ENHANCED: Now supports color vectors
  sort_within_groups = TRUE,
  annotation_title = "Gene Group"  # NEW: Custom title
)
```

#### `create_cell_annotations()`
```r
create_cell_annotations(
  exp_mat,                 # Required: Expression matrix
  percent_mat,             # Required: Percentage matrix
  split_pattern = "_",
  time_position = 1,
  celltype_start = 2,
  time_points_order = NULL,
  cell_types_order = NULL,
  time_color_palette = "Accent",    # ENHANCED: Now supports color vectors
  celltype_color_palette = "Dark2", # ENHANCED: Now supports color vectors
  show_time_annotation = TRUE,
  show_celltype_annotation = TRUE,
  time_point_title = "Time Point",  # NEW: Custom title
  cell_type_title = "Cell Type"     # NEW: Custom title
)
```

## Data Preparation

### Required Data Format

**Important**: This package expects your group identifiers to follow the format `"timepoint_celltype"` when you want to display both time and cell type information.

```r
library(dplyr)

# Method 1: Create combined group column (recommended for time course + cell type analysis)
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(group = paste(timepoint, celltype, sep = "_"))

# Example result: "Mock_Epidermis", "24hpi_Guard_cell", etc.

# Method 2: Use existing metadata columns for simpler analysis
# If you only want cell type information:
# group_by = "celltype"
# If you only want cluster information:
# group_by = "seurat_clusters"
```

## Usage Examples

### Basic Usage with All Default Parameters
```r
library(SingleCellComplexHeatMap)

# Minimal usage
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = gene_list
)
```

### Complete Parameter Customization
```r
# Full customization example
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = gene_list,
  gene_classification = gene_groups,
  group_by = "group",
  idents = c(9, 10, 13),
  time_points_order = c("Mock", "24hpi", "48hpi", "72hpi"),
  cell_types_order = c("Epidermis", "Guard_cell", "Mesophyll"),
  color_range = c(-2, 0, 3),
  color_palette = c("blue", "white", "red"),
  max_circle_size = 3,
  row_fontsize = 10,
  col_fontsize = 8,
  col_name_rotation = 45,
  row_title_fontsize = 12,
  col_title_fontsize = 12,
  show_heatmap_legend = TRUE,
  show_percentage_legend = TRUE,
  legend_side = "bottom",
  merge_legends = FALSE,
  percentage_legend_title = "% Expressing Cells",
  percentage_legend_labels = c("0", "20", "40", "60", "80", "100"),
  cell_border_color = "white",
  split_pattern = "_",
  gene_color_palette = "Dark2",
  time_color_palette = "Set3",
  celltype_color_palette = "Paired",
  show_gene_grouping = TRUE,
  show_time_annotation = TRUE,
  show_celltype_annotation = TRUE,
  split_by = "celltype",
  gene_group_title = "Gene Categories",
  time_point_title = "Treatment Time",
  cell_type_title = "Cell Types",
  show_cell_borders = TRUE,
  show_column_annotation = TRUE,
  gene_name_mapping = NULL
)
```

### Different Use Cases

#### 1. Full Analysis (Time + Cell Type + Gene Groups)
```r
# Prepare data with combined group
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(group = paste(sample, Cell_type, sep = "_"))

gene_groups <- list(
  "Senescence" = c("Gene1", "Gene2", "Gene3"),
  "Stress_Response" = c("Gene4", "Gene5", "Gene6")
)

heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"),
  gene_classification = gene_groups,
  group_by = "group",
  time_points_order = c("Mock", "24hpi", "48hpi", "72hpi"),
  cell_types_order = c("Cell_cycle", "Companion_cell", "Epidermis")
)
```

#### 2. Cell Type Only Analysis
```r
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = your_genes,
  group_by = "celltype",  # Use existing celltype column
  show_time_annotation = FALSE,
  split_by = "none"
)
```

#### 3. Simple Expression Analysis (No Grouping)
```r
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = your_genes,
  gene_classification = NULL,  # No gene grouping
  group_by = "seurat_clusters",
  show_time_annotation = FALSE,
  show_celltype_annotation = FALSE,
  split_by = "none"
)
```

#### 4. Custom Color Schemes
```r
# Using different color palettes
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = your_genes,
  gene_classification = gene_groups,
  color_range = c(-1.5, 0, 1.5, 3),  # 4-point color range
  color_palette = c("darkblue", "blue", "white", "red", "darkred"),
  gene_color_palette = "Spectral",
  time_color_palette = "Set2",
  celltype_color_palette = "Pastel1"
)
```

#### 5. Large Dataset Optimization
```r
# For large datasets, reduce visual complexity
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = your_genes,
  max_circle_size = 1.5,  # Smaller circles
  row_fontsize = 6,       # Smaller font
  col_fontsize = 6,
  cell_border_color = NA, # No borders for cleaner look
  merge_legends = TRUE    # Compact legends
)
```

#### 6. Publication-Ready Styling
```r
# High-quality publication figure
heatmap <- create_single_cell_complex_heatmap(
  seurat_object = seurat_obj,
  features = your_genes,
  gene_classification = gene_groups,
  color_range = c(-2, 0, 2),
  color_palette = c("#2166AC", "#F7F7F7", "#B2182B"),  # Custom colors
  max_circle_size = 2.5,
  row_fontsize = 10,
  col_fontsize = 10,
  row_title_fontsize = 12,
  col_title_fontsize = 12,
  percentage_legend_title = "Fraction of cells",
  percentage_legend_labels = c("0", "0.25", "0.50", "0.75", "1.0"),
  legend_side = "right",
  merge_legends = TRUE,
  gene_group_title = "Functional Categories",
  time_point_title = "Time Points",
  cell_type_title = "Cell Types",
  show_cell_borders = TRUE,
  cell_border_color = "white"
)
```

## Parameters

This section details all parameters for the `create_single_cell_complex_heatmap` function.

### Main Parameters
* `` `seurat_object` ``: A Seurat object containing single cell data. (No default)
* `` `features` ``: Character vector of gene names to plot. (No default)
* `` `gene_classification` ``: Named list where names are group labels and values are character vectors of gene names. (Default: `NULL`)
* `` `group_by` ``: Character string specifying the metadata column to group by. (Default: `"seurat_clusters"`)
* `` `idents` ``: Numeric or character vector specifying which cell groups to include. (Default: `NULL` for all)
* `` `time_points_order` ``: Character vector specifying order of time points. (Default: `NULL` for automatic ordering)
* `` `cell_types_order` ``: Character vector specifying order of cell types. (Default: `NULL` for automatic ordering)
* `` `color_range` ``: Numeric vector of length 3 or 4 specifying color mapping range for expression values. (Default: `c(-1, 0, 2)`)
* `` `color_palette` ``: Function or character vector specifying colors for the expression heatmap. (Default: `NULL`, then uses viridis palette if available, otherwise `color_palette_main`)
* `` `max_circle_size` ``: Numeric specifying maximum circle radius in mm for percentage representation. (Default: `2`)
* `` `row_fontsize` ``: Numeric specifying row name font size. (Default: `8`)
* `` `col_fontsize` ``: Numeric specifying column name font size. (Default: `9`)
* `` `col_name_rotation` ``: Numeric specifying column name rotation angle. (Default: `90`)
* `` `row_title_fontsize` ``: Numeric specifying row title font size (for gene groups). (Default: `10`)
* `` `col_title_fontsize` ``: Numeric specifying column title font size (for time/cell type splits). (Default: `10`)
* `` `show_heatmap_legend` ``: Logical indicating whether to show heatmap legend for expression values. (Default: `TRUE`)
* `` `show_percentage_legend` ``: Logical indicating whether to show percentage legend for circle sizes. (Default: `TRUE`)
* `` `legend_side` ``: Character string specifying legend position ("left", "right", "top", "bottom"). (Default: `"right"`)
* `` `cell_border_color` ``: Character string specifying cell border color. Use `NA` for no border. (Default: `"grey80"`)
* `` `split_pattern` ``: Character string used to split column names for parsing time points and cell types. (Default: `"_"` )
* `` `gene_color_palette` ``: **ENHANCED**: Character string specifying palette name OR character vector of colors for gene groups. (Default: `"Set1"`)
* `` `time_color_palette` ``: **ENHANCED**: Character string specifying palette name OR character vector of colors for time points. (Default: `"Accent"`)
* `` `celltype_color_palette` ``: **ENHANCED**: Character string specifying palette name OR character vector of colors for cell types. (Default: `"Dark2"`)
* `` `show_gene_grouping` ``: Logical indicating whether to show gene grouping annotation. (Default: `NULL`, resolves to `TRUE` if `gene_classification` is provided, else `FALSE`)
* `` `show_time_annotation` ``: Logical indicating whether to show time point annotation. (Default: `TRUE`)
* `` `show_celltype_annotation` ``: Logical indicating whether to show cell type annotation. (Default: `TRUE`)
* `` `split_by` ``: Character string specifying how to split columns: "time", "celltype", or "none". (Default: `"time"`)
* `` `merge_legends` ``: Logical indicating whether to merge legends if possible. (Default: `TRUE`)
* `` `percentage_legend_title` ``: Character string for percentage legend title. (Default: `"Expression %"`)
* `` `percentage_legend_labels` ``: Character vector for percentage legend labels. (Default: `c("0%", "25%", "50%", "75%", "100%")`)
* `` `return_data` ``: Logical indicating whether to return underlying data along with the heatmap object. (Default: `FALSE`)
* `` `save_plot` ``: Character string specifying file path to save plot (e.g., "my_heatmap.png"). (Default: `NULL`)
* `` `plot_width` ``: Numeric specifying plot width in inches when saving. (Default: `10`)
* `` `plot_height` ``: Numeric specifying plot height in inches when saving. (Default: `8`)
* `` `plot_dpi` ``: Numeric specifying plot resolution in DPI when saving. (Default: `300`)

### New Customization Parameters
* `` `gene_group_title` ``: **NEW**: Character string for gene group annotation title. (Default: `"Gene Group"`)
* `` `time_point_title` ``: **NEW**: Character string for time point annotation title. (Default: `"Time Point"`)
* `` `cell_type_title` ``: **NEW**: Character string for cell type annotation title. (Default: `"Cell Type"`)
* `` `show_cell_borders` ``: **NEW**: Logical indicating whether to show cell border lines. (Default: `TRUE`)
* `` `show_column_annotation` ``: **NEW**: Logical indicating whether to show column annotations. (Default: `TRUE`)
* `` `gene_name_mapping` ``: **NEW**: Named character vector for mapping gene names, where names are original gene names and values are display names. (Default: `NULL`)

### Data Source Control
* `` `assay` ``: Character string specifying which assay to use from Seurat object. (Default: `NULL`, uses active assay)
* `` `slot` ``: Character string specifying which slot to extract from assay (e.g., 'scale.data', 'data', 'counts'). (Default: `'scale.data'`)

### Clustering Control
* `` `cluster_cells` ``: Logical indicating whether to cluster cells/columns. (Default: `TRUE`)
* `` `cluster_features` ``: Logical indicating whether to cluster features/rows. (Default: `TRUE`)
* `` `clustering_distance_rows` ``: Character string specifying distance method for row clustering (e.g., "euclidean", "pearson"). (Default: `"euclidean"`)
* `` `clustering_distance_cols` ``: Character string specifying distance method for column clustering. (Default: `"euclidean"`)
* `` `clustering_method_rows` ``: Character string specifying clustering method for rows (e.g., "complete", "ward.D2"). (Default: `"complete"`)
* `` `clustering_method_cols` ``: Character string specifying clustering method for columns. (Default: `"complete"`)

### Color Mapping Control
* `` `color_palette_main` ``: Character vector of 3 colors for main heatmap gradient (used if `color_palette` is `NULL` and viridis is not available). (Default: `c("blue", "white", "red")`)
* `` `annotation_colors` ``: Named list specifying custom colors for annotations (e.g., `list(CellType = c(TypeA = "red"))`). (Default: `NULL`)

### Label and Legend Control
* `` `show_feature_names` ``: Logical indicating whether to show row names (gene names). (Default: `TRUE`)
* `` `feature_names_gp` ``: `gpar` object for controlling feature name appearance. (Default: `NULL`, then `gpar(fontsize = row_fontsize)`)
* `` `legend_title` ``: Character string for main heatmap legend title. (Default: `"Expression"`)

### Other Parameters
* `` `...` ``: Additional parameters passed to `ComplexHeatmap::Heatmap()` for maximum customization.

## Available Color Palettes

### RColorBrewer Palettes
- **Qualitative**: `"Set1"`, `"Set2"`, `"Set3"`, `"Pastel1"`, `"Pastel2"`, `"Paired"`, `"Dark2"`, `"Accent"`
- **Sequential**: `"Blues"`, `"Reds"`, `"Greens"`, `"Purples"`, `"Oranges"`, `"Greys"`
- **Diverging**: `"RdYlBu"`, `"RdBu"`, `"PiYG"`, `"PRGn"`, `"BrBG"`, `"PuOr"`, `"RdGy"`, `"Spectral"`

### Custom Color Examples
```r
# Viridis colors
color_palette = viridis::viridis(3)

# Custom RGB colors
color_palette = c("#440154", "#21908C", "#FDE725")

# Named colors
color_palette = c("navy", "white", "firebrick")
```

## Step-by-Step Workflow

```r
# Step 1: Prepare matrices
matrices <- prepare_expression_matrices(
  seurat_object = seurat_obj,
  features = your_genes,
  group_by = "group",
  idents = c(9, 10, 13)  # Specific cell clusters
)

# Step 2: Create gene annotations (if needed)
if (!is.null(gene_groups)) {
  gene_ann <- create_gene_annotations(
    exp_mat = matrices$exp_mat,
    percent_mat = matrices$percent_mat,
    gene_classification = gene_groups
  )
}

# Step 3: Create cell annotations
cell_ann <- create_cell_annotations(
  exp_mat = gene_ann$exp_mat_ordered,
  percent_mat = gene_ann$percent_mat_ordered,
  time_points_order = c("Mock", "24hpi", "48hpi", "72hpi"),
  cell_types_order = c("Epidermis", "Guard_cell", "Mesophyll")
)
```

## Column Naming Convention

For optimal functionality, ensure your column identifiers follow these patterns:

- **Time + Cell Type**: `"timepoint_celltype"` (e.g., "Mock_Epidermis", "24hpi_Guard_cell")
- **Complex naming**: `"timepoint_celltype_additional_info"` (e.g., "Mock_Guard_cell_cluster1")

The package automatically parses:
- Position 1: Time point
- Position 2+: Cell type (can include underscores)

## Requirements

- R (>= 4.0.0)
- ComplexHeatmap
- Seurat
- dplyr
- tidyr
- RColorBrewer
- circlize
- grid

## License

MIT License - see LICENSE file for details.
