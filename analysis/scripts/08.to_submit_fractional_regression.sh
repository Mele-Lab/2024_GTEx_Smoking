#!/bin/sh
#SBATCH --job-name=dsa
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/dsa_%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/dsa_%A_%a.err
#SBATCH --cpus-per-task=48
#SBATCH --array=41
#SBATCH --time=08:00:00
#SBATCH --constraint=highmem


module load gcc pcre2 R/4.2.0
#2-47
file=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/tissues_sorted.csv
export tissue=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file | cut -d "," -f2 | xargs)

echo $tissue 
Rscript /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/scripts/08.DSA.R -t $tissue


