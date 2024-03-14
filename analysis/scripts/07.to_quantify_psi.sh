#!/bin/sh
#SBATCH --job-name=PSI
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/err/%A.err
#SBATCH --cpus-per-task=20
#SBATCH --time=16:00:00

module purge
module load intel mkl python/3.10.2
path=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis

for tissue in $(awk -F ',' '{if(NR!=1) {print $2}}' ${path}/data/public/tissues_sorted.csv | tr -d \"); do 
	echo $tissue
	outpath=${path}/SUPPA/PSI_values/${tissue}/
	mkdir -p ${outpath}
	expression_file=${path}/SUPPA/TranscripExpressionFiles/${tissue}.transcript_TPM.txt
	for splicing_event in SE MX AF AL A5 A3 RI; do
		splicing_event_file=${path}/SUPPA/gencode.v26.annotation_events.ioe_${splicing_event}_strict.ioe
		output_file=${path}/SUPPA/PSI_values/${tissue}/${tissue}.${splicing_event}
		python3 /gpfs/projects/bsc83/utils/SUPPA-2.3/suppa.py psiPerEvent --ioe-file $splicing_event_file --expression-file $expression_file  -o $output_file --save_tpm_events 
	done
done





