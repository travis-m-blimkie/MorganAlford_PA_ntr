
# Libraries and data ------------------------------------------------------

library(tidyverse)


pa14_features <- read_tsv("../PA14_genome_info/pa14_operon_table.tsv") %>%
  select(Synonym, OperonID, Start, End, Strand) %>%
  rename(gene = "Synonym")




# Joining DE genes with features, operons, TSS ----------------------------

de_genes_1 <- de_genes %>%
  left_join(., pa14_features, by = c("gene" = "locus_tag"))

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
