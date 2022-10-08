# GDM-Project
Code summary of paper: Integrative Metagenomic  and Metabolomic Analysis Reveal Gut Microbiota-Derived Multiple Hits Connected to Development of Gestational Diabetes Mellitus in Humans



This is a R code record for GDM project
1. Differential_Abundance_Analysis.r is used for exploring the species with differentail abundance between GDM and NGT participants.

2. Diversity_Analysis.r is used for the exploration of beta diversity and alpha diversity.

3. PAPi_Score_Calculation.r is for calculation of PAPi Activity Score of metabolic pathways based on the metabolites abundance.

4. Bootstrap_Network_Analysis.r is for the bootstrap the samples in GDM and NGT groups. Then using the bootstrap tables, we constructed 100 subnetworks for two groups respectively and compared their network constructions.

5. RF_model_and_SHAP.r is used to build the RF model based on species relative abundance and apply SHAP to calculate the importance score of the species in the predicting the metabolic pathway.

6.  Selbal.r is used to explore the gut species balance between GDM and NGT.

7.  SPLS_Analysis.r aims to correlate the relative abundance of species, the relative abundance of metabolites, the PAPi score of metabolic pathways and the clinical phenotype.


If you find these codes are useful in your project, please consider to cite the paper.



Main packages used in this project are as follows:

Vegan	2.5-7	https://github.com/vegandevs/vegan/releases 

Selbal 0.1 https://github.com/malucalle/selbal 

Spls	2.2-3	https://cran.r-project.org/web/packages/spls/index.html 

PAPi	1.26.0 https://www.bioconductor.org/packages//2.12/bioc/html/PAPi.html 

Shapper	0.1.3 https://cran.r-project.org/web/packages/shapper/index.html 
