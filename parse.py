import os, sys
import re

key_files = []

with open('swissprot_psiblast.sam') as fil:
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

         target_dir = os.getcwd() + "/protein_fastas/" + key
         if(not os.path.exists(target_dir)):
            os.mkdir(target_dir)

         # save in target dir
         fil_path = target_dir + "/" + cur_line[0] + ".fasta"
         with open(fil_path, 'w') as new_file:
            new_file.write(">" + seq_name + "\n" + sequence)





