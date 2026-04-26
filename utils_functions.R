# ==============================================================================
# Script: Utility Functions
# Purpose: Common functions used across analyses
# ==============================================================================

# Function: Export Seurat object metadata
export_metadata <- function(seurat_obj, filename) {
  write.csv(seurat_obj@meta.data, file = filename, row.names = TRUE)
}

# Function: Quick QC plot
quick_qc_plot <- function(seurat_obj, output_file) {
  pdf(output_file, width = 12, height = 4)
  print(VlnPlot(seurat_obj, features = c("nFeature_Spatial", "nCount_Spatial", "mitoPercent"), 
                ncol = 3, pt.size = 0))
  dev.off()
}

# Function: Extract top markers
get_top_markers <- function(markers_df, n = 10) {
  markers_df %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
}

# Function: Save plots in multiple formats
save_plot_multi <- function(plot_obj, filename_prefix, width = 10, height = 8) {
  ggsave(paste0(filename_prefix, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(filename_prefix, ".png"), plot_obj, width = width, height = height, dpi = 300)
}

print("Utility functions loaded!")