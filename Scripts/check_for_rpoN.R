

# Load libraries ----------------------------------------------------------

library(tidyverse)


# Read in DE genes --------------------------------------------------------

de_genes <- list(
  lauren = read.csv("../DE_genes/de_lauren_pao1_edgeVScentre.csv", stringsAsFactors = F),
  shannon = read.csv("../DE_genes/de_shannon_pao1_swarmVswim_20190129.csv", stringsAsFactors = F)
)


# Read in rpoN genes ------------------------------------------------------

rpon_genes <- read.delim("../rpon_genes.tsv", stringsAsFactors = F, sep = "\t") %>%
  pull(locus_tag)


# Check if rpoN genes are DE ----------------------------------------------

de_rpon_genes <- de_genes %>%
  map(~filter(., locus_tag %in% rpon_genes))


# Save the DE rpoN genes --------------------------------------------------

map2(de_rpon_genes, names(de_rpon_genes), function(x, y)
  write_csv(x, path = paste0("../DE_genes/", y, "_pao1__de_rpoN_20190207.csv")))
