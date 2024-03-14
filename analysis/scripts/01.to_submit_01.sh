#!/bin/bash

#SBATCH --job-name=debug
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/%A.err
#SBATCH --qos=debug
#SBATCH --cpus-per-task=48

module load gcc pcre2 R/4.2.0

Rscript /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/scripts/01.Preparing_metadata.R


