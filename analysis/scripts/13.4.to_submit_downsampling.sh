#!/bin/bash

#SBATCH --job-name=downsampling
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/downsampling_%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/downsampling_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --array=1-4

module load R

export tissue=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Downsampling/tissues_methylation.txt | cut -f1)

echo $tissue
Rscript scripts/13.4.Model_downsampling.R -t $tissue




