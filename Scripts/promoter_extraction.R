
# Load required libraries -------------------------------------------------

library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(seqinr)
library(tidyverse)


# Read in the genome ------------------------------------------------------

pao1_genome <- readDNAStringSet(
  "../PAO1_genomic_info/pseudomonas_aeruginosa_PAO1_107.fna",
  format = "fasta"
  )


# Read in DE gene lists ---------------------------------------------------

master_lists <- list(
  lauren = read_csv("../Master_lists/masterList_lauren_20190129.csv"),
  shannon = read_csv("../Master_lists/masterList_shannon_20190129.csv")
) %>% map(~add_column(., name = rep("NC_002516")))


# Set up GRanges objects --------------------------------------------------

my_granges <- master_lists %>%
  map(~GRanges(
    .,
    seqnames = Rle(.$name),
    ranges = IRanges(start = .$promoter, end = .$position),
    strand = Rle(.$Strand_TSS),
    mcols = data.frame(
      locus_tag = .$locus_tag,
      operon = .$OperonID
    )
  ))


# Get upstream regions ----------------------------------------------------

my_promoters <- my_granges %>%
  map(~getSeq(pao1_genome, .) %>% as.data.frame())

my_promoters <- map2(master_lists, my_promoters, function(x, y)
  bind_cols(x, y)
  )


# Save the sequences ------------------------------------------------------

map2(my_promoters, names(my_promoters), function(a, b)
  write.fasta(as.list(a$x),
    names = a$locus_tag,
    file.out = paste0("../Promoter_sequences/promoterSeqs_", b, "_20190129.fna")
  ))
