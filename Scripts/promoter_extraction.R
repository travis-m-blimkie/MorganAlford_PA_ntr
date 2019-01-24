# Script to do promoter extraction from master gene lists


# Libraries and data ------------------------------------------------------

library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(seqinr)
library(tidyverse)

pao1_genome <- readDNAStringSet("../PAO1_genome/GCF_000006765.1_ASM676v1_genomic.fna", format = "fasta")


# DE genes ----------------------------------------------------------------

deg_master <- read.csv("Master_lists/DE_genes_master_list_20181001.csv")
deg_master$name <- rep("NC_002516.2")

deg_GR <- GRanges(
  seqnames = Rle(deg_master$name),
  ranges = IRanges(start = deg_master$promoter, end = deg_master$position),
  strand = Rle(deg_master$strand_TSS),
  mcols = data.frame(locus_tag = deg_master$gene,
                     operon = deg_master$OperonID)
)

deg_promoters <- getSeq(pao1_genome, deg_GR) %>% as.data.frame()

deg_promoters <- bind_cols(deg_master, deg_promoters)

write.fasta(as.list(deg_promoters$x),
            names = deg_promoters$gene,
            file.out = "../degenes_promoter_sequences_20181001.fasta")


# Gene clusters -----------------------------------------------------------

clusters_master <- read.csv("../clusters_master_list_20181001.csv")
clusters_master$name <- rep("NC_002516.2")

clusters_GR <- GRanges(
  seqnames = Rle(clusters_master$name), # Name must match header line of fasta file
  ranges = IRanges(start = clusters_master$promoter, end = clusters_master$position),
  strand = Rle(clusters_master$strand_TSS),
  mcols = data.frame(
    locus_tag = clusters_master$gene,
    operon = clusters_master$OperonID
  )
)

clusters_promoters <- getSeq(pao1_genome, clusters_GR) %>% as.data.frame()

clusters_promoters <- bind_cols(clusters_master, clusters_promoters)

write.fasta(as.list(clusters_promoters$x),
            names = clusters_promoters$gene,
            file.out = "../clusters_promoter_sequences_20181001.fasta")
