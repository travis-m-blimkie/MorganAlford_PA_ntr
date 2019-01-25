
# Load required libraries -------------------------------------------------

library(tidyverse)


# Read in PA14 gene info --------------------------------------------------

pa14_features <- read_tsv("../PA14_genome_info/pa14_operon_table.tsv") %>%
  select(Synonym, OperonID, Start, End, Strand) %>%
  rename(gene = "Synonym")


# Read in both DE gene lists ----------------------------------------------

de_genes <- list(
  lauren =  read_csv("../DE_genes/lauren_pa14_edgeVScentre.csv"),
  shannon = read_csv("../DE_genes/shannon_swarmVSswim_20190109.csv")
) %>% map(~select(., gene, log2FoldChange))


# Join DE genes and PA14 features -----------------------------------------

de_genes_info <- de_genes %>%
  map(~left_join(., pa14_features, by = "gene"))

de_genes_sorted <- de_genes_info %>%
  map(~arrange(., OperonID, Start))

de_genes_first <- de_genes_sorted %>%
  map(~distinct(., OperonID, .keep_all = T) %>% drop_na(Start))


# Add promoter location ---------------------------------------------------

de_genes_master <- de_genes_first %>%
  map(~mutate(., Promoter = Start - 250) %>%
        select(gene, OperonID, Promoter, Start, End, Strand, log2FoldChange))


# Save the results --------------------------------------------------------

map2(de_genes_master, names(de_genes_master), function(x, y)
  write_csv(x, path = paste0("../DE_genes/masterList_", y, "_20190125.csv")))
