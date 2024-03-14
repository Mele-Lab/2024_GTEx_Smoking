#!/bin/bash

#SBATCH --job-name=tpm
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/%A.err
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00


module load gcc pcre2 R/4.2.0

Rscript /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/scripts/07.Event_coordinates.R



