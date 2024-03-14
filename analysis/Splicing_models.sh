#Extracting tpm values per tissue
Rscript scripts/06.Preparing_splicing_metadata.R #If HPC: sbatch scripts/06.to_submit_splicing_metadata.sh

#Extract PSI values using SUPPA based on annotation. See: https://github.com/comprna/SUPPA#installation. We generate annotaton files "ioe" for each type of splicing event
module load intel mkl python/3.10.2 #We load the modules we use in our HPC system to run python3
python3 /gpfs/projects/bsc83/utils/SUPPA-2.3/suppa.py generateEvents -i data/public/gencode.v26.annotation.gtf -o SUPPA/gencode.v26.annotation_events.ioe -f ioe -e SE SS MX RI FL

sbatch scripts/07.to_quantify_psi.sh #Get PSI values per event type

#Merging PSI values per tissue
ls tissues/ | while read tissue; 
do scripts/07.Creating_PSI_files.R -t $tissue
done

sbatch scripts/07.to_submit_event_coordinate.sh

#Running Differential Splicing analysis model + hierarchical partitioning
sbatch scripts/08.to_submit_fractional_regression.sh


#Analyzing changes in functional domains:
Rscript 08.Parsing_for_Pfam.R

#git clone https://github.com/Mele-Lab/2022_GTExTranscriptome_CellGenomics_fromSplicingEventsToProteinDomains.git
sbatch run-Prot-Dom-Mapper-SLURM.sh #with edits in nextflow.config

Rscript scripts/08.differential_splicing_analysis.R #Getting output object file to share
Rscript scripts/08.differential_splicing_numbers.R #Analyzing results
