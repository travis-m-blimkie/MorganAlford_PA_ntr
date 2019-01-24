
# Libraries and data ------------------------------------------------------

library(tidyverse)

pao1_features <- read_tsv("PAO1_genome/GCF_000006765.1_ASM676v1_feature_table.txt") %>%
  filter(feature == "gene") %>%
  rename(., length = feature_interval_length) %>%
  select(locus_tag, start, end, strand, length)

pao1_operons <- read_tsv("PAO1_genome/pao1_operons_20180703.opr") %>%
  select(Synonym, OperonID, Product)

pao1_tss <- read_csv("PAO1_genome/allTSS_forlab.csv") %>%
  select(locus_tag, TSS, position, distance)

de_genes <- read_csv("Original_lists/DE_relAspoTvsPAO1_planktonic_original.csv") %>%
  select(gene, FC)


# Joining DE genes with features, operons, TSS ----------------------------

de_genes_1 <- de_genes %>%
  left_join(., pao1_features, by = c("gene" = "locus_tag")) %>%
  left_join(., pao1_operons, by = c("gene" = "Synonym")) %>%
  left_join(., pao1_tss, by = c("gene" = "locus_tag"))

de_genes_2 <- de_genes_1 %>%
  arrange(OperonID, start, position)

de_genes_3 <- de_genes_2[!duplicated(de_genes_2$OperonID), ]

de_genes_4 <- de_genes_3 %>% drop_na(start)


#' Adds start coordinate to "position" column from TSS info, if that cell is NA
#' This way "start" information is contained in a single column for all genes,
#' including those without specific TSS information

for (i in 1:nrow(de_genes_4)) {
  if (is.na(de_genes_4[i, "position"])) {
    de_genes_4[i, "position"] = de_genes_4[i, "start"]
  }
}


# Changing strand based on TSS type ---------------------------------------

de_genes_5 <- de_genes_4
de_genes_5$strand_TSS <- de_genes_5$strand

for (i in 1:nrow(de_genes_5)) {
  if (de_genes_5$TSS[i] %in% c("antisense", "primary antisense")) {
    if (de_genes_5$strand[i] == "+") {
      de_genes_5$strand_TSS[i] <- "-"
    } else if (de_genes_5$strand[i] == "-") {
      de_genes_5$strand_TSS[i] <- "+"
    }
  }
}


# Selecting columns of interest, saving the result ------------------------

de_genes_master <- de_genes_5 %>%
  mutate(promoter = position - 250) %>%
  select(gene, promoter, position, start, end, length, strand, strand_TSS, TSS, OperonID, Product, FC)

write_csv(de_genes_master, path = "../DE_genes_master_list_20181001.csv")
