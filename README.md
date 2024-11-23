# NCYdiversity

My scripts and data for the NCYdiversity (Non-Conventional Yeast diversity) project.

* bioconductor_pip_ITS_functions.R functions used by the modified bioconductor pipeline for ITS. 

* bioconductor_pip_v7_4_0224_ITS.R the pipeline itself mostly taken from https://benjjneb.github.io/dada2/ITS_workflow.html but with several modifications. 

* primer_pairs_fungi.txt a service table containing most frequently used primer pairs for ITS amplification, work in progress  

* bigdata folder contains three scripts modeled after https://benjjneb.github.io/dada2/bigdata.html, with several modifications; you should 

  +  run the getHeader script first (to partition the sequences in groups of manageable size)
  
  + run the bigdata_part1 script as many times as the groups you created in the sequences with getHeader, to do preprocessing, filtering, estimation of the error model, ASV inference and merging and produce sequence tables. 
  
  + run the bigdata_part2 script to merge the sequence tables, do taxonomy assignment, infer phylogenetic tree, assemble a phyloseq object and produce files for imponrt in FoodMicrobionet  
  
* FMBN5_0_1 a folder containing R lists for FoodMicrobionet, table specifications and summary statistics

## Acknowledgements.  

This work was carried out within the PRIN 2022 PNRR Project NCY diversity P20229JMMH and received funding from the European Union Next-GenerationEU, CUP C53D23007560001 (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2,  INVESTIMENTO 1.4 – D.D. 1048 14/07/2023). This script and its contents reflects only the authors’ views and opinions,  neither the European Union nor the European Commission can be considered  responsible for them.
