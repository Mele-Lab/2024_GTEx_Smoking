#!/bin/bash


mkdir -p output/DEGS_enrichement/BP.up
mkdir -p output/DEGS_enrichement/BP.down




orsum.py --gmt  data/gmt/hsapiens.GO_BP.name.gmt --files 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/AdiposeSubcutaneous_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/EsophagusMucosa_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/Lung_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/Pancreas_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/SkinSunExposedLowerleg_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/Stomach_BP_up.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.upregulated/Thyroid_BP_up.csv' --fileAliases  AdiposeSubcutaneous EsophagusMucosa Lung Pancreas SkinSunExposedLowerleg Stomach Thyroid --outputFolder output/DEGS_enrichement/BP.up --minTermSize 1
orsum.py --gmt  data/gmt/hsapiens.GO_BP.name.gmt --files 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/ArteryTibial_BP_down.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/EsophagusMucosa_BP_down.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/Lung_BP_down.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/Pancreas_BP_down.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/Stomach_BP_down.csv' 'output/DEGS_enrichement/enriched_terms.go.bp.downregulated/Thyroid_BP_down.csv' --fileAliases  ArteryTibial EsophagusMucosa Lung Pancreas Stomach Thyroid --outputFolder output/DEGS_enrichement/BP.down --minTermSize 1
orsum.py -v > orsum.version
