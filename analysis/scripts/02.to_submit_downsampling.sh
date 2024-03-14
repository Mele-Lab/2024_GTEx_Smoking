#!/bin/bash

#SBATCH --job-name=downsampling
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/downsampling_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/downsampling_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --qos=debug

module load R

cat Downsampling/tissues_4.txt | while read tissue; 
do 
echo $tissue
for i in {1..50}; do
echo $i
scripts/02.voom_limma.Smoking.R -t $tissue -d -m /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/Downsampling/11/$tissue/metadata_downsampling_$i.rds
scripts/02.voom_limma.Smoking.R -t $tissue -d -m /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/Downsampling/25/$tissue/metadata_downsampling_$i.rds
scripts/02.voom_limma.Smoking.R -t $tissue -d -m /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/Downsampling/35/$tissue/metadata_downsampling_$i.rds
done; done



