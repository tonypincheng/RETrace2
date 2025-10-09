#!/usr/bin/env Rscript

################################################################################
# Circular Phylogenetic Tree Visualization with Outer Rings
#
# This script generates circular phylogenetic trees with continuous outer ring
# annotations for the RETrace2 MSH2 mouse data. It creates two ring layers:
# - Inner ring: Tissue/group colors
# - Outer ring: Cell type colors
#
# Features:
# - Circular layout with customizable tip labels (plain or colored backgrounds)
# - Two continuous color rings for hierarchical annotations
#
# Usage:
#   Rscript view_phylo_ggtree.R \
#     --samplesheet <samples.csv> \
#     --tree_file <tree.nwk> \
#     --prefix <output_prefix> \
#     [--color_background]
#
# Input requirements:
# - Newick tree file
# - CSV samplesheet with columns: sample_id, color, group, 
#   cell_type_color, cell_type_assignment
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
  library(ggtree)
  library(ggtreeExtra)
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggnewscale)
})

# Argument parser
parser <- ArgumentParser(description = 'Visualize phylogenetic tree using ggtree with outer ring')
parser$add_argument('--samplesheet', required = TRUE, help = 'Path to sample information CSV file')
parser$add_argument('--tree_file', required = TRUE, help = 'Path to newick tree file')
parser$add_argument('--prefix', required = TRUE, help = 'Output prefix for visualization files')
parser$add_argument('--color_background', action = 'store_true', default = FALSE, help = 'Apply colors to text background instead of node circles')

args <- parser$parse_args()

# Load tree
tree <- read.tree(args$tree_file)

# Load sample data
samples <- read.csv(args$samplesheet, stringsAsFactors = FALSE)

# Handle missing values with defaults
if (!'color' %in% colnames(samples)) {
  samples$color <- 'grey'
} else {
  samples$color <- ifelse(is.na(samples$color) | samples$color == '', 'grey', samples$color)
}

if (!'cell_type_color' %in% colnames(samples)) {
  samples$cell_type_color <- 'grey'
} else {
  samples$cell_type_color <- ifelse(is.na(samples$cell_type_color) | samples$cell_type_color == '', 'grey', samples$cell_type_color)
}

# Prepare data for tree - only color column for tip labels
tip_label_data <- samples %>% 
  select(sample_id, color) %>%
  rename(label = sample_id)
rownames(tip_label_data) <- tip_label_data$label

# Prepare data for inner ring with group names for legend
group_ring_data <- samples %>%
  select(sample_id, color, group) %>%
  rename(label = sample_id)

# Prepare data for outer ring with cell type names for legend
cell_type_ring_data <- samples %>%
  select(sample_id, cell_type_color, cell_type_assignment) %>%
  rename(label = sample_id)

# Custom rename dictionary for group names
rename_dict <- list(
  'Brain_LeftAnteriorCortex' = 'Brain Left Anterior Cortex',
  'Brain_LeftPosteriorCortex' = 'Brain Left Posterior Cortex',
  'Kidney_Left-Cortex-1' = 'Kidney Left Cortex',
  'Kidney_Right-Cortex-1' = 'Kidney Right Cortex',
  'Liver_LeftLobe' = 'Liver Left Lobe',
  'Liver_MedianLobe' = 'Liver Median Lobe',
  'Liver_RightLobe' = 'Liver Right Lobe'
)

# Function to rename based on dictionary, leave unchanged if not in dict
rename_labels <- function(x, dict) {
  sapply(x, function(name) {
    if (name %in% names(dict)) {
      dict[[name]]
    } else {
      name
    }
  })
}

# Apply custom renaming to group names
group_names_clean <- rename_labels(unique(samples$group), rename_dict)
group_colors_map <- setNames(unique(samples$color), group_names_clean)

# Update the group_ring_data to use cleaned names
group_ring_data$group <- rename_labels(group_ring_data$group, rename_dict)

# Cell type names - replace underscores with spaces
cell_type_names_clean <- gsub("_", " ", unique(samples$cell_type_assignment))
cell_type_colors_map <- setNames(unique(samples$cell_type_color), cell_type_names_clean)

# Update the cell_type_ring_data to use cleaned names
cell_type_ring_data$cell_type_assignment <- gsub("_", " ", cell_type_ring_data$cell_type_assignment)

# Basic tree plot - circular layout only
p <- ggtree(tree, layout = 'circular', branch.length = 'none', size = 1)

# Join ONLY color data to tree for tip labels
p <- p %<+% tip_label_data

# Add tip labels with or without background color
if (args$color_background) {
  # Text with colored background - smaller font to avoid overlap
  p <- p + geom_tiplab(aes(fill = color), 
                       size = 2.5, color = 'black', 
                       geom = 'label', 
                       label.padding = unit(0.1, "lines")) +
    scale_fill_identity()
} else {
  # Colored node circles with plain text - smaller font with offset
  p <- p + geom_tippoint(aes(color = color), size = 3) +
    geom_tiplab(size = 2.5, offset = 0.3) +  # Small offset to push text away from nodes
    scale_color_identity()
}

# Add inner continuous ring for group color with legend
p <- p + new_scale_fill() + 
  geom_fruit(
    data = group_ring_data,
    geom = geom_tile,
    mapping = aes(y = label, x = 2, fill = group),
    pwidth = 0.1,
    width = 3,
    offset = 0.72
  ) + 
  scale_fill_manual(
    values = group_colors_map,
    name = "Tissue (Inner Ring)",
    guide = guide_legend(order = 1)  # Order 1 = top
  )

# Add outer continuous ring for cell_type_color with legend
p <- p + new_scale_fill() + 
  geom_fruit(
    data = cell_type_ring_data,
    geom = geom_tile,
    mapping = aes(y = label, x = 2, fill = cell_type_assignment),
    pwidth = 0.001,
    width = 3,
    offset = 0.3
  ) + 
  scale_fill_manual(
    values = cell_type_colors_map,
    name = "Cell Type (Outer Ring)",
    guide = guide_legend(order = 2)  # Order 2 = bottom
  )

# Adjust plot margins and theme, and position legends
# Minimize white space with minimal margins
p <- p + theme(
  plot.margin = margin(0, 0, 0, 0),  # Remove all margins
  legend.box = "vertical",
  legend.text = element_text(size = 20),  # Larger legend text
  legend.title = element_text(size = 21, face = "bold"),  # Larger, bold legend titles
  legend.background = element_rect(fill = "white", colour = NA)  # White background, no border
) +
  theme(
    legend.position.inside = c(0.85, 0.5),  # Apply position separately
    legend.justification = c(0.5, 0.5)  # Center anchor
  )

# Determine output size - wider than tall to accommodate circular tree + legend
width <- 18
height <- 14

# Save outputs
ggsave(paste0(args$prefix, '.viewPhylo.pdf'), p, width = width, height = height, dpi = 500)
ggsave(paste0(args$prefix, '.viewPhylo.png'), p, width = width, height = height, dpi = 500)

cat("Tree visualization saved to:\n")
cat(paste0("  ", args$prefix, ".viewPhylo.pdf\n"))
cat(paste0("  ", args$prefix, ".viewPhylo.png\n"))

