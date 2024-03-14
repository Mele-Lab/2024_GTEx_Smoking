#Preparing metadata
Rscript scripts/01.Preparing_metadata.R

#Running models for each clinical trait per tissue to update the metadata (to keep or exclude the clinical traits for downstream analysis)
ls tissues/ | while read tissue; 
do scripts/02.voom_limma.clinical_traits.R -t $tissue
done

#Get Tabl S1 with the diseases included per tissue
Rscript scripts/02.Table_S1.R

#Running models with smoking per tissue
ls tissues/ | while read tissue; 
do scripts/02.voom_limma.Smoking.R -t $tissue
done

#We submit the hierarchical partitioning to our SLURM HPC machine, but it is not mandatory, we can directly call Rscript /scripts/03.Hier_part.R
sbatch 03.to_submit_hier_part.sh

#Code to get zenodo tables from the differential gene expression analysis results 
Rscript scripts/04.Parsing_expression_results.R #We also have the option to run sbatch 04.to_submit_parsing.sh if SLURM is avalable

#Running interactions
ls tissues/ | while read tissue; 
do scripts/02.voom_limma.Smoking.R -t $tissue -i "Age:Smoking BMI:Smoking Sex:Smoking Ancestry:Smoking"
done

#Preparing code for Figure 2
Rscript 04.Demographic_traits.R

#Downsampling
Rscript 02.Downsampling.R #To create metada
sbatch scripts/02.to_submit_downsampling.sh #Create metadata and models
sbatch scripts/02.to_submit_downsampling_2.sh #Compute medians
