
# Load required libraries -------------------------------------------------

library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(seqinr)
library(tidyverse)


# Read in the genome ------------------------------------------------------

pa14_genome <- readDNAStringSet(
  "../PA14_genome_info/Pseudomonas_aeruginosa_UCBPP-PA14_109.fna",
  format = "fasta"
  )


# Read in DE gene lists ---------------------------------------------------

master_lists <- list(
  lauren = read_csv("../DE_genes/masterList_lauren_20190125.csv"),
  shannon = read_csv("../DE_genes/masterList_shannon_20190125.csv")
) %>% map(~add_column(., name = rep("NC_008463")))


# Set up GRanges objects --------------------------------------------------

my_granges <- master_lists %>%
  map(~GRanges(
    .,
    seqnames = Rle(.$name),
    ranges = IRanges(start = .$Promoter, end = .$Start),
    strand = Rle(.$Strand),
    mcols = data.frame(
      locus_tag = .$gene,
      operon = .$OperonID
    )
  ))


# Get upstream regions ----------------------------------------------------

my_promoters <- my_granges %>%
  map(~getSeq(pa14_genome, .) %>% as.data.frame())

my_promoters <- map2(master_lists, my_promoters, function(x, y)
  bind_cols(x, y)
  )


# Save the sequences ------------------------------------------------------

map2(my_promoters, names(my_promoters), function(a, b)
  write.fasta(as.list(a$x),
    names = a$gene,
    file.out = paste0("promoterSeqs_", b, "_20190125.fasta")
  ))
