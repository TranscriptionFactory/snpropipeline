# Get gene names from genomic coordinates of SNPs associated interactions at 
# protein-protein interfaces


# Input: <Protein Code of Interactor> <chromosome>:<coordinates>:<Ref>:<Var>

setwd("/ix/djishnu/Aaron")

library(parallel)
library(foreach)
library(iterators)
library(doParallel)

p.cluster = parallel::makeCluster(parallel::detectCores() - 2)

doParallel::registerDoParallel(p.cluster)

library(tidyverse)
library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# load ppi SNPs
ppi_input = read.delim(paste0(getwd(), "/data/rsids.hg38.e-3.HuHu.interface.txt"),
                       header = F, col.names = c("protein", "coords"))

# separate columns of txt file into columns in dataframe
ppi_input = separate(ppi_input, col = coords, into = c("chr", "coord", "ref", "var"), sep = ":")

# convert coordinates to list in form chr:start-end (start and end are the same)
ppi_input$search_coords = paste0(ppi_input$chr, ":", ppi_input$coord, "-", ppi_input$coord)

# Get overlapping exons (in parallel) - Retyr s ~17k results
convertExonsToGenes_Parallel = function(ppi_input_row) {
  # Function to map SNP to exons by using GRanges.
  # GRanges will find exons that overlap with the SNP (only want protein coding) and will return them
  
  # Returns annotations of ensembl gene ID, entrez ID, gene name and sequencing strand
  # Using in parallel, each core needs a copy of the EnsDb (see below)
  exons_found = exonsByOverlaps(EnsDb.Hsapiens.v86, ranges = GRanges(ppi_input_row$search_coords),
                                filter = list(TxBiotypeFilter(list("protein_coding",
                                                                   "processed_transcript"))),
                                columns = c("gene_id", "entrezid", "gene_name",
                                            "seq_strand"))
  
  genes_found = exons_found@elementMetadata
  
  genes_found$strand = exons_found@strand@values
  
  ppi_input_row = cbind(ppi_input_row, genes_found)
  
  return(ppi_input_row)
}

# (dont use) Get overlapping transcripts (in parallel) - returns ~40k results
getTranscripts_from_SNPs_Parallel = function(ppi_input_row) {
  # Function to map SNP to exons by using GRanges.
  # GRanges will find transcripts  that overlap with the SNP (only want protein coding) and will return them
  
  # Returns annotations of ensembl gene ID, entrez ID, gene name and sequencing strand
  # Using in parallel, each core needs a copy of the EnsDb (see below)
  tx_found = ensembldb::transcriptsByOverlaps(EnsDb.Hsapiens.v86, ranges = GRanges(ppi_input_row$search_coords),
                                # columns = c("gene_id", "tx_id"))
                                filter = list(TxBiotypeFilter(list("protein_coding",
                                                                   "processed_transcript"))),
                                columns = c("gene_id", "entrezid", "gene_name"))
  
  genes_found = tx_found@elementMetadata
  ppi_input_row = cbind(ppi_input_row, genes_found)
  
  return(ppi_input_row)
}

# prepare workers for parallel computing
# each core needs a copy of all of the libraries and the EnsDb
# (every core starts with a blank R environment so we need to initialize
# each environment and export any variables we want accessed)
clusterExport(p.cluster, c("ppi_input"), envir = environment())
clusterEvalQ(p.cluster, library("ensembldb"))
clusterEvalQ(p.cluster, library("tidyverse"))
clusterEvalQ(p.cluster, library("biomaRt"))
clusterEvalQ(p.cluster, library("EnsDb.Hsapiens.v86"))
clusterEvalQ(p.cluster, library("GenomicFeatures"))

# start cluster (pass EnsDb again here...clusterEvalQ calls prepare the individual
# cores and the above
# the .packages calls are for the forreach package
results = foreach(i = 1:nrow(ppi_input),
                               .combine = rbind,
                               .packages = c("ensembldb", "EnsDb.Hsapiens.v86", "GenomicFeatures")) %dopar% {

                                 convertExonsToGenes_Parallel(ppi_input[i,])
                               }

parallel::stopCluster(p.cluster)

results = results %>% dplyr::distinct()

results.distinct = data.frame(results) %>% dplyr::filter(entrezid != "NA")

results.distinct = results.distinct[, names(results.distinct) != "exon_id"]

results.distinct = distinct(results.distinct)

unique_genes = data.frame(unique(results.distinct$gene_name))

# write results
# write_tsv(unique_genes, file = paste0(getwd(), "/data/gene_list.tsv"))
# write_tsv(results, file = paste0(getwd(), "/data/gene_annotations_full.tsv"))

