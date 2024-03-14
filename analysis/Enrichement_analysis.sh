#Run enrichement analysis (this also generates the orsum script)
Rscript scripts/05_enrichement_analysis.R

# Summarise the terms using ORSUM
sh scripts/05_run_orsum.sh

# Generate the data for the figure
Rscript scrips/enrichement_results.R