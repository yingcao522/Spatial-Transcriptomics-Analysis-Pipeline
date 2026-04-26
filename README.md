# Spatial-Transcriptomics-Analysis-Pipeline
This repository contains a comprehensive analysis pipeline for spatial transcriptomics data integration with single-cell RNA-seq.

## Structure
- `00_setup_environment.R` - Environment setup
- `01_data_loading_QC.R` - Data loading and quality control
- `02_normalization_integration.R` - SCTransform and Harmony integration
- `03_dimensionality_reduction_clustering.R` - UMAP/tSNE and clustering
- `04_cell_type_annotation.R` - Cell type annotation
- `05_spatial_scRNA_integration.R` - **Spatial deconvolution (Anchor-based)**
- `06_differential_expression_analysis.R` - DEG analysis
- `07_functional_enrichment.R` - GO/KEGG enrichment
- `08_cellchat_analysis.R` - Cell-cell communication
- `09_spatial_domain_analysis.R` - Spatial region analysis
- `10_gene_set_scoring.R` - AUCell gene set scoring
- `11_visualization.R` - Publication figures
- `12_statistical_analysis.R` - Statistical tests
- `utils_functions.R` - Utility functions

## Usage
Run scripts sequentially from 00 to 12. Adjust parameters in each script according to your data.

## Requirements
See `00_setup_environment.R` for required packages.
