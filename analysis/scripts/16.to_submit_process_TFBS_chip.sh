#!/bin/bash

#SBATCH --job-name=tfbs
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/tfbs_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/tfbs_%A.err
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00

module load R

Rscript scripts/16.Process_TFBS_ChipAtlas.R




