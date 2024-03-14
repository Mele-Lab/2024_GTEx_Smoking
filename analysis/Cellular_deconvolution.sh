# Prepare data 
Rscript scripts/09_preProcess_scRNA_seq.R

# Run deconvolution 
Rscript scripts/09_single_cell_deconvolution.R

# Compare macrophage populations
Rscript scripts/09_compare_macrophages.R