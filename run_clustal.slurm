#!/bin/bash
#SBATCH --job-name=ppi_clustal
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aar126@pitt.edu
#SBATCH --qos=short
#SBATCH --time=0-23:01:00


#before running, run this in the terminal one time to create a file in each folder
# with all the sequences in one .fasta (run it from the /data/fastas_clustal/protein_fastas folder
# or deal with writing the directory names
# to replace the files:
# delete first
# for file in *; do rm ${file}/clustal_in.fasta; done

#for file in *; do for f in ${file}/*.fasta; do (cat "${f}"; echo) >> ${file}/clustal_in.fasta; done; done
module load python/ondemand-jupyter-python3.8

source activate clustal_env

file_dir=/ix/djishnu/Aaron/data/fastas_clustal/protein_fastas

for file in ${file_dir}/*; do
  # if ! [[ -f ${file}_50ITER_clustal_output ]]; then
    clustalo -i ${file}/clustal_in.fasta --dealign -o ${file}_mBed_50_clustal_output --iterations=50 --outfmt=fa;
  # fi
  # if ! [[ -f ${file}_50ITER_clustal_output ]]; then
  #   clustalo -i ${file}/clustal_in.fasta --dealign -o ${file}_50ITER_clustal_output --full --iterations=25;
  # fi
done
    