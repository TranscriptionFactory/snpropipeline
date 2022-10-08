library(bio3d) # read.fasta
library(tidyverse)

setwd("/ix/djishnu/Aaron/scripts/snpropipeline/")
source("entropy_func.R")

# already calculated
if(F){
  
  output_folder = file.path("/ix/djishnu/Aaron/data/fastas_clustal/clustal_output_50_iter/")
  
  output_folder_files = list.files(output_folder)
  
  entropy_list = list()
  alignment_list = list()
  # try with first one
  for(file in output_folder_files) {
    fsta = bio3d::read.fasta(paste0(output_folder, file))
    en_val = entropy(fsta)
    entropy_list = append(entropy_list, list(en_val))
    alignment_list = append(alignment_list, list(fsta$ali))
    
    # entropy_list = append(entropy_list, list("entropy" = entropy(fsta), "alignments" = fsta$ali))
  }
  
  
  
  saveRDS(entropy_list, "/ix/djishnu/Aaron/data/output_data/entropy_values.RDS")
  saveRDS(alignment_list, "/ix/djishnu/Aaron/data/output_data/alignment_entropy_values.RDS")
}

entropy_list = readRDS("/ix/djishnu/Aaron/data/output_data/entropy_values.RDS")

mut_index = readRDS("/ix/djishnu/Aaron/data/output_data/mutant_sequence_index.RDS")

df = data.frame()

for(i in 1:length(alignment_list)){
  index = i
  seqnames = row.names(alignment_list[[i]])
  goi = grep('^ENSP', seqnames)
  if(length(goi > 0)){
    for(j in 1:length(goi)){
      parsed = separate(data.frame(seqnames[goi]), col = 1, into = c("ensembl_protein_id", "ensembl_transcript_id", "chr", "coord", "id"))
      parsed$index = i
      df = rbind(df, parsed)
    }
  }
}
df = df %>% drop_na()

df$chr = as.double(df$chr)
df$coord = as.double(df$coord)
# merge with residue start #s
df = df %>% left_join(mut_index, by = c("chr", "coord"))

H_annotated_SNPs = df %>% select(-id) %>% distinct()
# get H values at position in index of entropy file
H_annotated_SNPs$H = -1
for(i in 1:nrow(H_annotated_SNPs)){
  ind = H_annotated_SNPs[i,]$index
  Hvals = entropy_list[[ind]]$H
  
  
  H_annotated_SNPs[i,"H"] = Hvals[H_annotated_SNPs[i,]$residue_start]
}

saveRDS(H_annotated_SNPs, "/ix/djishnu/Aaron/data/output_data/ENTROPY_VALS_FOR_SNPS.RDS")

# parse sequence names: fasta files are named with tx#_protein#_chr_pos_ind
# where index is ithe index in the mutations df
#look up in this table