library(tidyverse)
library(rentrez)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(ensembldb)

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#load ensembl annotations
setwd("/ix/djishnu/Aaron/data")

ens_tbl = read_tsv("/ix/djishnu/Aaron/data/gene_annotations_full.tsv")

# get unique ensembl gene IDs

ens_genes = ens_tbl %>% distinct(gene_id)

# get transcript (primary tx)
filt = listFilters(ensembl)[54, 1]

attr = listAttributes(ensembl)[c(18, 32,3, 1, 22), 1]

tx = getBM(attributes = attr, filters = filt, values = ens_genes$gene_id, mart = ensembl)

# all mane select txs
tx_select = tx %>% dplyr::filter(transcript_mane_select != '')

# for filtering by transcript biotype (just use mane select transcripts instead)
# tx_bt = unique(tx$transcript_tsl)
# tx_bt = tx_bt[grepl('tsl1', tx_bt)]
#
#filter for only tsl1 and protein coding
# txs = tx %>% dplyr::filter(transcript_biotype == "protein_coding", transcript_tsl %in% tx_bt)
# 
# ens_distinct = txs %>%
#   left_join(ens_genes, by = c("ensembl_gene_id" = "gene_id"))

# map transcripts to sequences

# get transcript sequences
tx_seq = getBM(attributes = c(attr[3],"coding", "ccds"),
               filters = listFilters(ensembl)[56,1], 
               values = tx_select$ensembl_transcript_id, mart = ensembl)


tx_seq = tx_seq %>% left_join(tx_select, by = "ensembl_transcript_id")

# left join with large gene annotation table

# remove exons (we don't need duplicated rows because of different exon annotations)
tx_seqs_by_pos = (ens_tbl %>% dplyr::select(-entrezid, -tx_biotype, -strand, -exon_id)) %>% distinct()
tx_seqs_by_pos = tx_seqs_by_pos %>% left_join(tx_seq, by = c("gene_id" = "ensembl_gene_id"))

# remove na columns because they aren't annotated with high quality transcripts
tx_seqs_by_pos = tx_seqs_by_pos %>% drop_na()

# we can also remove the MANE select version so that we keep dimensions small

tx_seqs_by_pos = separate(tx_seqs_by_pos, col = "transcript_mane_select", into = c("transcript_mane_select", "mane_version"),
              sep = "\\.") %>% dplyr::select(-mane_version)

saveRDS(tx_seqs_by_pos, file = "/ix/djishnu/Aaron/data/output/MANE_select_tx.RDS")
tx_seq = readRDS("/ix/djishnu/Aaron/data/output/MANE_select_tx.RDS")

#map coordinate to transcript position

#coordinate ranges
# coord_range = GRanges(paste0(data$chr, ":", data$coord, "-", (data$coord+1) ))
coord_range = GRanges(paste0(tx_seq$chr, ":", tx_seq$coord))

# map genome coordinates to transcript coordinates (this takes a while to run)
gnm_tx = genomeToTranscript(coord_range, EnsDb.Hsapiens.v86)

# saveRDS(gnm_tx, "/ix/djishnu/Aaron/data/output/transcript_ranges.RDS")
gnm_tx = readRDS("/ix/djishnu/Aaron/data/output/transcript_ranges.RDS")

#unlist GRanges object and convert to dataframe
# main ones need to be removed are exon id/rank because there are multiple mappings of
# different exons to same transcript, so a join creates many, many rows (~2 million)
tx_ranges = as.data.frame(gnm_tx@unlistData) %>% 
  dplyr::select(-names, -seq_strand, -exon_id, -exon_rank)

# join by chromosomal position
tx_ranges$search_coords = paste0(tx_ranges$seq_name, ":", tx_ranges$seq_start, "-", tx_ranges$seq_end)

# join with transcripts found earlier (already selected for high quality)
unique_tx_select = unique(tx_seq$transcript_mane_select)
#get ensembl transcripts associated with these
unique_tx_select_ensembl = (tx_seq %>% dplyr::filter(transcript_mane_select %in% unique_tx_select))$ensembl_transcript_id

# filter returned transcripts
tx_ranges_filtered = tx_ranges %>% dplyr::filter(tx_id %in% unique_tx_select_ensembl)

# these are all the variants (excluding the specific ppis)
txs_joined_ranges = (tx_seq %>% dplyr::select(-protein)) %>% 
  left_join(tx_ranges_filtered, by = c("ensembl_transcript_id" = "tx_id", "search_coords" = "search_coords"),
            keep = F) %>%
  dplyr::filter(start != '') %>% distinct() %>% drop_na() %>% dplyr::select(-seq_start, -seq_end, -seq_name)

# saveRDS(txs_joined_ranges, "/ix/djishnu/Aaron/data/output/mapped_transcript_sequences.RDS")
txs_joined_ranges = readRDS("/ix/djishnu/Aaron/data/mapped_transcript_sequences.RDS")


# get protein sequence for transcript
# ranges_uncomp = IRanges(pos = txs_joined_ranges$start, 
#                         # width = txs_joined_ranges$width, 
#                         names = txs_joined_ranges$names)
ranges_uncomp = IRanges(start = txs_joined_ranges$start,
                        width = txs_joined_ranges$width,
                        names = txs_joined_ranges$ensembl_transcript_id)

cds_pos_from_transcript = transcriptToCds(ranges_uncomp, db = EnsDb.Hsapiens.v86)
# saveRDS(cds_pos_from_transcript,"/ix/djishnu/Aaron/data/output/cds_pos_from_transcript.RDS")

# cds_pos_from_transcript = readRDS("/ix/djishnu/Aaron/data/output/cds_pos_from_transcript.RDS")

# map transcript to protein 
# cds_ranges = IRanges(start = cds_pos_from_transcript@start,
#                      width = cds_pos_from_transcript@width,
#                      names = cds_pos_from_transcript@NAMES)
# 

protein_from_tx = transcriptToProtein(ranges_uncomp, db = EnsDb.Hsapiens.v86)

protein_from_txdf = as.data.frame(protein_from_tx) %>% distinct() %>% dplyr::filter(cds_ok == "TRUE")
# saveRDS(protein_from_tx, "/ix/djishnu/Aaron/data/output/protein_transcript_map.RDS")

names(protein_from_txdf)[1:2] = c("residue_start", "residue_end")
names(protein_from_txdf)[4] = "ensembl_protein_id"

# join everything

fully_joined_ids = (txs_joined_ranges %>% dplyr::select(-width)) %>% left_join(protein_from_txdf, 
                                                   by = c("ensembl_transcript_id" = "tx_id",
                                                   "start" = "tx_start",
                                                   "end" = "tx_end"), keep = F)
# saveRDS(fully_joined_ids, "/ix/djishnu/Aaron/data/output/mapped_snp_to_protein.RDS")
