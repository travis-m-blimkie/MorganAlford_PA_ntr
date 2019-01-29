
# Load required libraries -------------------------------------------------

library(tidyverse)


# Read in PAO1 gene info --------------------------------------------------

pao1_features <- read_tsv("../PAO1_genomic_info/pseudomonas_aeruginosa_PAO1_operons.txt") %>%
  select(Synonym, OperonID, Start, End, Strand) %>%
  rename(locus_tag = "Synonym")

pao1_tss <- read_csv("../PAO1_genomic_info/erinGill_pao1_tss.csv") %>%
  select(locus_tag, TSS, position, distance)


# Read in both DE gene lists ----------------------------------------------

de_genes <- list(
  lauren =  read_csv("../DE_genes/de_lauren_pao1.csv"),
  shannon = read_csv("../DE_genes/de_shannon_pao1.csv")
) %>% map(~select(., locus_tag, log2FoldChange))


# Join DE genes and PA14 features -----------------------------------------

de_genes_joined <- de_genes %>%
  map(~plyr::join_all(list(., pao1_features, pao1_tss), by = "locus_tag", type = "left") %>%
        arrange(., OperonID, Start, position)
  )


de_genes_first <- de_genes_joined %>%
  map(~distinct(., OperonID, .keep_all = T) %>% drop_na(., ... = Start))


# Getting TSS info --------------------------------------------------------

#' Adds start coordinate to "position" column from TSS info, if that cell is NA
#' This way "start" information is contained in a single column for all genes,
#' including those without specific TSS information
add_pos <- function(x) {
  for (i in 1:nrow(x)) {
    if (is.na(x[i, "position"])) {
      x[i, "position"] <- x[i, "Start"]
    }
  }
  return(x)
}

de_genes_pos <- de_genes_first %>% map(~add_pos(.))


# Changing strand based on TSS type ---------------------------------------

de_genes_tss <- de_genes_pos %>% map(~mutate(., Strand_TSS = Strand))

# Function to change strand based on TSS type from Gill et al
change_strand <- function(x) {
  for (i in 1:nrow(x)) {
    if (x$TSS[i] %in% c("antisense", "primary antisense")) {
      if (x$Strand[i] == "+") {
        x$Strand_TSS[i] <- "-"
      } else if (x$Strand[i] == "-") {
        x$Strand_TSS[i] <- "+"
      }
    }
  }
  return(x)
}

de_genes_tss <- de_genes_tss %>% map(~change_strand(.))


# Selecting columns of interest -------------------------------------------

de_genes_master <- de_genes_tss %>%
  map(~mutate(., promoter = position - 250) %>%
        select(locus_tag, promoter, position, Start, End, Strand, Strand_TSS, TSS, OperonID)
  )


# Save the result ---------------------------------------------------------

# map2(de_genes_master, names(de_genes_master), function(x, y)
#   write_csv(x, path = paste0("../Master_lists/masterList_", y, "_20190129.csv")))
