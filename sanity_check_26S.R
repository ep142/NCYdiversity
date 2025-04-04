# sanity check SILVA LSU

# the script checks if the DADA2 trainsets created from the SIVA v138.2 SSU
# reference database
# SILVA_138.2_LSURef_NR99_tax_silva.fasta downloaded from
# https://www.arb-silva.de/no_cache/download/archive/release_138_2/Exports/

# the sequences used for testing were downloaded from:
# https://migale.pages.mia.inra.fr/metabarfood/benchmark.html#Composition_of_mock_communities


library(phylotools)
library(tidyverse)
library(tidysq)
library(dada2)
library(caret)

# read the .csv file
DNA_mock <- read_csv(file.path("seq_mock_migale","METABARFOOD - Results on mock communities - Choice of the bioinformatic solution.csv"))
# remove the sequences with non standard nucleotides
DNA_mock <- DNA_mock  |> 
  mutate(gsubbed = str_remove_all(D1.D2, "[^ACGTacgt]")) |>
  mutate(seq_identical = (gsubbed == D1.D2)) |>
  dplyr::filter(seq_identical)

seqs <- DNA_mock$D1.D2
names(seqs) <- DNA_mock$Strain

# assign taxonomy
ref_fasta <- "SILVA_LSUfungi_nr99_v138_2_toGenus_trainset.fasta"
assign_species_ref <- "SILVA_LSUfungi_assignSpecies.fasta"
taxtab <- assignTaxonomy(seqs, refFasta = ref_fasta, multithread = TRUE)
taxtab2 <- addSpecies(taxtab, assign_species_ref) 
taxtab2 <- as.data.frame(taxtab2) |>
  rownames_to_column(var = "D1.D2") |>
  mutate(species = if_else(!is.na(Species), str_c(Genus, Species, sep = " "), Species))

refseq_26S_assgn <- bind_cols(DNA_mock,
                              select(taxtab2, genus_assgnd = Genus, species_assgnd = species))
xtabs(~Genus + genus_assgnd, data = refseq_26S_assgn)

matches_long <- refseq_26S_assgn |>
  group_by(Genus, genus_assgnd) %>%
  tally() 
# matching ide (non very accurate)
matches_count <- matches_long |>
  mutate(matched = str_detect(genus_assgnd, Genus)) |>
  mutate(matched = if_else(is.na(matched), F, matched))
matched_TF <- matches_count |>
  ungroup() |>
  summarise(assignment = sum(n), .by = c(Genus, matched))

matches_wide <- matched_TF |>
  pivot_wider(id_cols = Genus, names_from = matched, 
              values_from = assignment, values_fill = 0)


# let's try with confusion matrix
# this is to handle the use od clades in assignment
refseq_26S_assgn <- refseq_26S_assgn |> 
  mutate(genus_assgnd_2 = if_else(str_detect(genus_assgnd, Genus), Genus, genus_assgnd)) |>
  mutate(genus_assgnd_2 = if_else(is.na(genus_assgnd_2), "No match", genus_assgnd_2)) |>
  mutate(genus_assgnd_2 = str_remove_all(genus_assgnd_2, "[\\[\\]]"))

genus_levels <- sort(unique(c(refseq_26S_assgn$Genus, refseq_26S_assgn$genus_assgnd_2)))

# a new df
refseq_26S_assgn_2 <- refseq_26S_assgn |>
  select(truth = Genus, predicted = genus_assgnd_2) |>
  mutate(truth = factor(truth, levels = genus_levels),
         predicted = factor(predicted, levels = genus_levels))


conf_matrix_LSU <- confusionMatrix(refseq_26S_assgn_2$truth,refseq_26S_assgn_2$predicted)
conf_matrix_LSU

results_by_class <- as.data.frame(conf_matrix_LSU$byClass)

sink(file = "confusionMatrix_LSU,txt")
confusionMatrix(refseq_26S_assgn_2$truth,refseq_26S_assgn_2$predicted)
sink(file = NULL)

# Acknowledgements  
# This work was was carried out within the PRIN 2022 PNRR Project NCY diversity P20229JMMH and received funding from the European Union Next-GenerationEU, CUP C53D23007560001 (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2,  INVESTIMENTO 1.4 – D.D. 1048 14/07/2023). This script and its contents reflects only the authors’ views and opinions,  neither the European Union nor the European Commission can be considered  responsible for them.