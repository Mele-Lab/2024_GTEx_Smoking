#!/bin/bash

#SBATCH --job-name=downsampling
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/comp_downsampling_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/comp_downsampling_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00


module load R

Rscript scripts/02.Downsampling_2.R



