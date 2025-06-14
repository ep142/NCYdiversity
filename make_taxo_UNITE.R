# UNITE genus
# a script for creating lineages from genus up from UNITE trainsets
require(phylotools)
require(tidyverse)

# open the taxonomic reference
taxo_db_fn <- file.choose()
# "sh_general_release_dynamic_s_04.04.2024.fasta"

taxo_df <- read.fasta(taxo_db_fn)
# need to split the taxonomy
taxo_UNITE <- taxo_df |>
  separate_wider_delim(cols = seq.name, names = c("accession","lineage"),
                       cols_remove = F, delim = "k__") |>
  separate_wider_delim(cols = lineage, names = c("Kingdom","lineage"),
                       delim = ";p__") |>
  separate_wider_delim(cols = lineage, names = c("Phylum","lineage"),
                       delim = ";c__") |>
  separate_wider_delim(cols = lineage, names = c("Class","lineage"),
                       delim = ";o__") |>
  separate_wider_delim(cols = lineage, names = c("Order","lineage"),
                       delim = ";f__") |>
  separate_wider_delim(cols = lineage, names = c("Family","lineage"),
                       delim = ";g__") |>
  separate_wider_delim(cols = lineage, names = c("Genus","Species"),
                       delim = ";s__")
write_tsv(taxo_UNITE, file = "taxo_UNITE_gen_rel_04_04_2024.txt")
# unique genera

UNITE_genera <- taxo_UNITE |>
  select(Kingdom:Genus) |>
  distinct() |>
  arrange(Genus)
any(duplicated(UNITE_genera$Genus))
# duplicated_genera <- UNITE_genera$Genus[which(duplicated(UNITE_genera$Genus))]

write_tsv(UNITE_genera, file = "UNITE_genus_lineage.txt")
