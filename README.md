# NYCdiversity

My scripts and data for the NYCdiversity project.

* bioconductor_pip_ITS_functions.R functions used by the modified bioconductor pipeline for ITS. 

* bioconductor_pip_v7_4_0224_ITS.R the pipeline itself mostly taken from https://benjjneb.github.io/dada2/ITS_workflow.html but with several modifications. 

* primer_pairs_fungi.txt a service table containing most frequently used primer pairs for ITS amplification, work in progress  

* bigdata folder contains three scripts modeled after https://benjjneb.github.io/dada2/bigdata.html, with several modifications; you should 

  +  run the getHeader script first (to partition the sequences in groups of manageable size)
  
  + run the bigdata_part1 script as many times as the groups you created in the sequences with getHeader, to do preprocessing, filtering, estimation of the error model, ASV inference and merging and produce sequence tables. 
  
  + run the bigdata_part2 script to merge the sequence tables, do taxonomy assignment, infer phylogenetic tree, assemble a phyloseq object and produce files for imponrt in FoodMicrobionet
