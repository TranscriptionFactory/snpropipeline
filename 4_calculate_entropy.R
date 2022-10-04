library(bio3d) # read.fasta
library(tidyverse)

source("entropy_func.R")

output_folder = file.path("/ix/djishnu/Aaron/data/fastas_clustal/clustal_output/")

output_folder_files = list.files(output_folder)

# try with first one
f = read.fasta(paste0(output_folder, output_folder_files[1]))

entropy(f)
