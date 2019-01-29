

# Load libraries ----------------------------------------------------------

library(tidyverse)


# Read in ortholog info ---------------------------------------------------

orthos <- read_tsv("../PAO1_genomic_info/pa14_pao1_orthologs.txt")


# Read in DE gene lists ---------------------------------------------------

de_genes_pa14 <- list(
  lauren = read_csv("../DE_genes/de_lauren_pa14_edgeVcentre.csv"),
  shannon = read_csv("../DE_genes/de_shannon_pa14_swarmVswim_20190109.csv")
)


# Left join to map orthologs ----------------------------------------------

de_genes_pao1 <- de_genes_pa14 %>%
  map(~left_join(., orthos, by = c("gene" = "Query_Gene")) %>%
        drop_na(Hit_Gene) %>%
        select(Hit_Gene, 2:9) %>%
        rename(locus_tag = "Hit_Gene")
  )


# Save the new PAO1 genes -------------------------------------------------

map2(de_genes_pao1, names(de_genes_pao1), function(x, y)
  write_csv(x, path = paste0("../DE_genes/de_", y, "_pao1.csv")))
