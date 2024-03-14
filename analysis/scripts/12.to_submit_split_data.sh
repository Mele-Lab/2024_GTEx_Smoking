#!/bin/bash

#SBATCH --job-name=split_methylation
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/%A.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --constraint=highmem

module load gcc pcre2 R/4.2.0

Rscript scripts/12.Split_data.R
