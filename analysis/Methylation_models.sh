sbatch scripts/12.to_submit_split_data.sh

Rscript 12.PEER_correlations.R

for tissue in BreastMammaryTissue ColonTransverse KidneyCortex Lung MuscleSkeletal Ovary Prostate Testis WholeBlood
do
sbatch scripts/12.to_submit_model.sh $tissue
done

Rscript scripts/13.1.Methylation_results.R #6A

sbatch scripts/13.1.Replication_to_literature_plot.R

#Studying enrichment of DMPs in particular genomics locations
Rscript scripts/13.2.location_of_DMPs_ChromHMM.R
Rscript scripts/13.2.location_of_DMPs_EPIC.R

#Studying functional enrichments of DMPs
Rscript scripts/13.2.functional_enrichments_DMP.R

#Exploring methylation and expression correlations
13.3.to_submit_correlations.sh #which calls "Rscript scripts/13.3.Correlation.R"

#Exploring overlap between smoking and age
Rscript scripts/14.Methylation_age.R

#Exploring reversible effects on DNA methylation after smoking cessation
Rscript scripts/14.Methylation_reversibility.R


#Downsampling to model the exact same samples in both expression and methylation
Rscript scripts/15.Subset_metadata.R

for tissue in BreastMammaryTissue ColonTransverse Lung MuscleSkeletal Ovary Prostate Testis WholeBlood
do
Rscript scripts/02.voom_limma.Smoking.R -t $tissue -s
Rscript scripts/12.Model.R -t $tissue -s
done

Rscript scripts/15.Subset_reversibility.R

#Downsampling and permuting:
Rscript 13.4.Methylation_downsampling.R #Create metadata
sbatch scripts/13.4.to_submit_downsampling.sh #Run models

#TFBS
Rscript 16.Process_TFBS_ChipAtlas.R
Rscript 16.TFBS_ChipAtlas_github.R

#Other analyses
Rscript 17.general_methylation.R
Rscript 17.bivalent_enhancers_correlations.R

