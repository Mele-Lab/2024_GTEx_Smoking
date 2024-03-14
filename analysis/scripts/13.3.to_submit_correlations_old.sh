#!/bin/bash

#SBATCH --job-name=cor
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/cor_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/cor_%A.err
#SBATCH --cpus-per-task=48
#SBATCH --time=08:00:00

module load gcc pcre2 R/4.2.0

Rscript scripts/13.3.Correlation_methylation_expression.R
