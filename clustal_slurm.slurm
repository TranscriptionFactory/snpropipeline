#!/bin/bash

#SBATCH --job-name=ppi_clustal
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --qos=normal


#before running, run this in the terminal one time to create a file in each folder
# with all the sequences in one .fasta (run it from the /data/fastas_clustal/protein_fastas folder
# or deal with writing the directory names

#for file in *; do for f in ${file}/*.fasta; do (cat "${f}"; echo) >> ${file}/clustal_in.fasta; done; done
module load python/ondemand-jupyter-python3.8

source activate clustal_env

file_dir=/ix/djishnu/Aaron/data/fastas_clustal/protein_fastas

for file in ${file_dir}/*; do
clustalo -i ${file}/clustal_in.fasta --dealign -o ${file}_clustal_output --full --iterations=5;
done
    