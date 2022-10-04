import os, sys
import re
import argparse

parser = argparse.ArgumentParser(description='Parse SAM files into FASTA')
parser.add_argument("--blast",required=True,help="SAM file from BLAST")
parser.add_argument("--path",required=True,help="Path to save FASTA files (separated by gene id)")

args = parser.parse_args()

key_files = []

with open(str(args.blast)) as fil:
   lines = fil.readlines()

   for index in range(0, len(lines)):
      if str(lines[index]).startswith('@'):
         # split string
         cur_line = re.split(r'(\t)+|(SN:)', str(lines[index]))
         key_files.append(str(cur_line[6]))
      else:
         cur_line = str(lines[index]).replace(".", "_") # remove .
         cur_line = cur_line.split() #whitespace
         # match with key file
         key = cur_line[2]
         seq_name = cur_line[0] + "_" + key

         # the entire sequence is in this index
         sequence = cur_line[9]

         target_dir = str(args.path) + key
         if(not os.path.exists(target_dir)):
            os.mkdir(target_dir)

         # save in target dir
         fil_path = target_dir + "/" + cur_line[0] + ".fasta"
         with open(fil_path, 'w') as new_file:
            new_file.write(">" + seq_name + "\n" + sequence)





