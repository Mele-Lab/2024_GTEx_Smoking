#!/bin/bash

#SBATCH --job-name=DML
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/meth_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/meth_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --qos=debug

module load gcc pcre2 R/4.2.0

echo $1
Rscript /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/scripts/12.Model.R -t $1 -s
