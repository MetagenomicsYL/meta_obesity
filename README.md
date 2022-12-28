# meta_obesity
Meta analysis of Obesity

## work flow
You can follow the instructions below to get the main results in the manuscript:

- There are several commonly used functions in the **function_library.R**. This file should be imported during the differential abundance analysis for taxa and function modules.
 
- The phylum-level composition in each study was generated using **phyla_upset_plot.R**

- **alpha_div_metrics.R** was used to calculate alpha diversity metrics and F/B ratio between obese and healthy individuals in different studies.

- **combined_alpha_diversity.R** was used to calculate the combined alpha diversity using the random effect model and generate the forest plot of alpha diversity.

- **beta_diversity.R** was used for the PCoA analysis for obesity.

- The Stratified Wilcoxon test for the microbial taxa and functional pathways was calculated according to **preprocess_taxa_table.R** and **process_function_table.R**

- The random effect model for each genus was calculated using **select_microbial_genera.R** and the volcano plot for the identification of core taxa were generated using **volcano_plot_REM_Wilcoxon.R**

- **diversity_BMI_ethnicity.R**, **diversity_BMI_AGP_HMP.R**, and **species_BMI_correlation.R** were used to compute the correlation between BMI with alpha diversity and microbial taxa.

- **pathway_abundance_analysis.R** was used to identify the differentially abundant functional modules between obese and healthy controls and correlation_function_genera.R was used to correlate the taxa and functional modules.
