library(bio3d) # read.fasta
library(tidyverse)

source("entropy_func.R")

output_folder = file.path("/ix/djishnu/Aaron/data/fastas_clustal/clustal_output/")

output_folder_files = list.files(output_folder)

entropy_list = list()
# try with first one
for(file in output_folder_files) {
  fsta = read.fasta(paste0(output_folder, file))
  entropy_list = append(entropy_list, list("entropy" = entropy(fsta), "alignments" = fsta$ali))
}

saveRDS(entropy_list, "/ix/djishnu/Aaron/data/output_data/entropy_values.RDS")
