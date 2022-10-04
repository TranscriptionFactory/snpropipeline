library(tidyverse)
library(rentrez)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(Biostrings)
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#load ensembl annotations
setwd("/ix/djishnu/Aaron/data")

data = readRDS("/ix/djishnu/Aaron/data/output/mapped_snp_to_protein.RDS") %>% drop_na()

names(data)[which(names(data) == "start")] = "tx_start"
names(data)[which(names(data) == "end")] = "tx_end"
if(T){
#### joining cds with tx coordinates for mapping
cds_pos = readRDS("/ix/djishnu/Aaron/data/output/cds_pos_from_transcript.RDS") %>% as.data.frame() %>% distinct()
names(cds_pos)[1:2] = c("cds_start", "cds_end")

data = (data %>% dplyr::select(-width)) %>% left_join(cds_pos, by = c("ensembl_transcript_id" = "names",
                                          "tx_start" = "tx_start",
                                          "tx_end" = "tx_end"))

#### getting peptide seqs
unique_proteins = unique(data$ensembl_protein_id)
# data = data %>% dplyr::filter( (coord == seq_start),(start == tx_start),(end == tx_end))
# map these for seqs
protein_seqs = getBM(attributes = c("peptide", "ensembl_peptide_id"), filters = "ensembl_peptide_id",
                     values = unique_proteins, mart = ensembl)

data = data %>% left_join(protein_seqs, by = c("ensembl_protein_id" = "ensembl_peptide_id"))
}

data = separate(data, col = "peptide", sep = "\\*", into = c("peptide", "asterisk")) %>% dplyr::select(-asterisk)

substitute_base = function(rw){
  tx = Biostrings::DNAString(rw$coding)
  # make edit
  if(tx[rw$tx_start] == rw$ref){
    tx[rw$tx_start] = rw$var
  
    rw$mut_tx = toString(tx)
  }
  return(rw)
}


index = data %>% dplyr::select(coding, cds_start, ref, var, 
                               ensembl_protein_id, peptide, gene_id, ensembl_transcript_id)

mut_txs = sapply(index$coding, Biostrings::DNAString)
mut_prots = sapply(index$peptide, Biostrings::AAString)
for(i in 1:length(mut_txs)){
  sq = mut_txs[[i]]
  pos = index[i,]$cds_start
  letter = DNAString(index[i,]$var)
  if(length(letter) > 1){
    inds = seq(pos, pos+length(letter) - 1, by = 1)
    replaceLetterAt(sq, at = inds, letter = letter)
  } else {
    replaceLetterAt(sq, at = pos, letter = letter)
    
  }
  mut_txs[[i]] = sq
  mut_prots[[i]] = translate(sq, if.fuzzy.codon = "solve")
}

names(mut_prots) = index$ensembl_protein_id
# create mutated consensus
# unique_proteins = unique(data$ensembl_protein_id)
# 
# data_consensus_proteins = data
# for(i in 1:length(unique_proteins)){
#   selected_proteins = data_consensus_proteins %>% filter(ensembl_protein_id == unique_proteins[i])
#   
# }


# saveRDS(index, "/ix/djishnu/Aaron/data/index_file.RDS")
# saveRDS(mut_txs, "/ix/djishnu/Aaron/data/mutant_transcripts.RDS")    
# saveRDS(mut_prots, "/ix/djishnu/Aaron/data/output/mutant_proteins.RDS")
# index = readRDS("/ix/djishnu/Aaron/data/output/index_file.RDS")
# 
# index_single_proteins = (data %>% dplyr::select(gene_id, ensembl_transcript_id, ensembl_protein_id, peptide)) %>% distinct()
# 
# prots_string_set = AAStringSet(index_single_proteins$peptide)
# names(prots_string_set) = index_single_proteins$ensembl_protein_id
# 
# saveRDS(single_proteins, "/ix/djishnu/Aaron/data/output/WT_peptide_list.RDS")
# writeXStringSet(prots_string_set, filepath = "/ix/djishnu/Aaron/data/output/protein_list.fasta", format = "fasta")

# proteins = Biostrings::readAAStringSet("/ix/djishnu/Aaron/data/output/protein_list.fasta")

# save each set of mutated proteins (by gene) into folder for later parsing
protein_names = unique(index$ensembl_protein_id)

mut_prots_AAset = AAStringSet(mut_prots)

mut_prots_df = data.frame(mut_prots_AAset@ranges)

mut_prots_df$peptide = unlist(lapply(mut_prots_AAset, Biostrings::toString))
mut_prots_df$tx_id = index$ensembl_transcript_id

# remove asterisk, but doesn't matter - will be parsed in python
mut_prots_df = separate(mut_prots_df, col = "peptide", sep = "\\*", into = c("peptide", "asterisk")) %>% dplyr::select(-asterisk)
# saveRDS(mut_prots_df, "/ix/djishnu/Aaron/data/output/mutant_protein_sequences.RDS")


for(i in 1:length(protein_names)) {
  # get proteins
  same_protein = mut_prots_df %>% dplyr::filter(names == protein_names[i])
  
  # write these to folder
  protein_dir = paste0("/ix/djishnu/Aaron/data/output/fastas_by_protein/", protein_names[i])
  
  if(!file.exists(protein_dir)){
    dir.create(protein_dir)
  }
  protein_string_set = AAStringSet(same_protein$peptide)
  names(protein_string_set) = paste0(same_protein$names, "_", same_protein$tx_id)
  
  writeXStringSet(protein_string_set, paste0(protein_dir, "/", protein_names[i],".fasta"), format = "fasta")
  ref_pep = AAStringSet(index$peptide[i])
  names(ref_pep) = paste0(protein_names[i], "_ref")
  writeXStringSet(ref_pep, paste0(protein_dir, "/", protein_names[i],"_ref.fasta"), format = "fasta")
}












