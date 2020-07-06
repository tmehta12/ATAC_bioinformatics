#!/usr/bin/env python3

# Here's a Python 3 script that splits each full MEME file into multiple files based on the each motif delimiter. Script ran like:
# python3 /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/split_meme.py -i input_file -o output_dir -t file_type -s species

# Example input file:

# MEME version 4
#
# ALPHABET= ACGT
#
# strands: + -
#
# Background letter frequencies (from uniform background):
# A 0.25000 C 0.25000 G 0.25000 T 0.25000
#
# MOTIF TP53
#
# letter-probability matrix: alength= 4 w= 18 nsites= 18978 E= 0
#   0.433264        0.059557        0.376924        0.130255
#   0.761994        0.008407        0.176330        0.053269
#   0.002186        0.956520        0.000632        0.040663
#   0.882014        0.015081        0.048274        0.054631
#   0.080989        0.017713        0.006489        0.894809
#   0.000392        0.000098        0.999510        0.000000
#   0.012530        0.616957        0.005597        0.364916
#   0.056923        0.801279        0.019330        0.122468
#   0.060438        0.626814        0.113539        0.199208
#   0.221360        0.084727        0.572168        0.121745
#   0.151020        0.031856        0.751980        0.065144
#   0.440121        0.001611        0.546496        0.011771
#   0.000000        0.999742        0.000207        0.000052
#   0.892179        0.017257        0.010837        0.079727
#   0.091727        0.005931        0.004473        0.897869
#   0.015589        0.001573        0.979597        0.003242
#   0.039609        0.150415        0.005227        0.804749
#   0.071456        0.445777        0.040049        0.442718
#
# MOTIF RELA
#
# letter-probability matrix: alength= 4 w= 10 nsites= 18 E= 0
#   0.000000        0.222222        0.611111        0.166667
#   0.000000        0.000000        0.944444        0.055556
#   0.000000        0.000000        1.000000        0.000000
#   0.611111        0.000000        0.388889        0.000000
#   0.555556        0.166667        0.222222        0.055556
#   0.111111        0.000000        0.000000        0.888889
#   0.000000        0.000000        0.000000        1.000000
#   0.000000        0.111111        0.000000        0.888889
#   0.000000        1.000000        0.000000        0.000000
#   0.000000        1.000000        0.000000        0.000000


import os # get paths
import argparse # provide parameter based system arguments
import time # this is used for time delays

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True, help="input filename")
parser.add_argument("-o", "--output_dir", required=True, help="output directory")
parser.add_argument("-t", "--file_type", required=True, help="the file origin used for outfile naming e.g. CS or CW or JP")
parser.add_argument("-s", "--species", required=True, help="species output ID e.g. mz, pn, ab, nb, on")

args = parser.parse_args()

# read the input file
with open(args.input_file, 'r') as input_file1:
    input_data = input_file1.read()

# with open("2c_JASPAR_ab.meme", 'r') as input_file1:
#     input_data = input_file1.read()

# store the MEME header as head
header = ""
with open(args.input_file, 'r') as input_file1:
    for i in range(8):
        header += input_file1.readline()

# header = ""
# with open("2c_JASPAR_ab.meme", 'r') as input_file1:
#     for i in range(8):
#         header += input_file1.readline()

with open(args.input_file, 'r') as input_file2: # open the input file
    input_data = input_file2.read() # read each line
    for i,part in enumerate(input_data.split("MOTIF")): # split each section by the string motif and assign the line number 'i' with the rest of the information as 'part'
        # print(part.split(None, 1)[0]) # this will use no delimiter to get the very first word in each part string i.e. the motif ID - use this for naming
        output_filename = part.split(None, 1)[0] + "_" + args.file_type + "_" + args.species + ".meme"
        output_path = os.path.join(args.output_dir, output_filename)
        with open(output_path, 'w') as output_file: # create a newfile for each section (i)
            output_file.write(header + '\n' + "MOTIF" + part) # write the MEME header, the section and

# with open("2c_JASPAR_ab.meme", 'r') as input_file2:
#     input_data = input_file2.read()
#     for i,part in enumerate(input_data.split("MOTIF")):
#         # print(part.split(None, 1)[0]) # this will use no delimiter to get the very first word in each part string i.e. the motif ID - use this for naming
#         with open(part.split(None, 1)[0] + "_" + "JP" + ".meme", 'w') as newfile:
#             newfile.write(header + '\n' + "MOTIF" + part)

