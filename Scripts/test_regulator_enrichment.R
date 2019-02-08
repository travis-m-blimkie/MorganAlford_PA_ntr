
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(tRavis)


# Read in DE genes --------------------------------------------------------

de_genes <- list(
  lauren = read.csv("../DE_genes/de_lauren_pao1_edgeVScentre.csv", stringsAsFactors = F),
  shannon = read.csv("../DE_genes/de_shannon_pao1_swarmVswim_20190129.csv", stringsAsFactors = F)
) %>% map(~pull(., locus_tag))


# Read in list of regulators ----------------------------------------------

regulators_df <- bind_rows(read.csv("../all_regulon_genes_unique.csv", stringsAsFactors = F),
                           read.delim("../rpon_genes.tsv", sep = "\t", stringsAsFactors = F)
)

regulators_ls <- split(regulators_df$locus_tag, regulators_df$TF)


# Test for enrichment -----------------------------------------------------

enrich_raw <- map(de_genes, function(a)
  map(regulators_ls, function(y) tr_test_enrichment(query_genes = a, enrichment_set = y, total_genes = 5688))
  )


# Convert result to a data frame ------------------------------------------

enrich_df <- enrich_raw %>%
  map(~tibble(
    Regulator = names(.),
    PValue = as.numeric(.[names(.)])
  ))


# Correct for multiple testing and filter ---------------------------------

enrich_filter <- enrich_df %>%
  map(~mutate(., PAdj = p.adjust(PValue, method = "BH")) %>%
        filter(PAdj <= 0.05)
     )

enrich_one_df <- bind_rows(enrich_filter, .id = "Source")


# Save the combined output ------------------------------------------------

write_csv(enrich_one_df, path = "../regulator_enrichment_result_20190207.csv")
