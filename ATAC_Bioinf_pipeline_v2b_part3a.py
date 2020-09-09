#!/usr/bin/python

################################################################################################################
## ---------------------------
##
## Script name: ATAC_Bioinf_pipeline_v2b_part3a.R
##
## Purpose of script: Filtering blastn output for potential hits to mitochondrial genome
##
## Author: Dr. Tarang K. Mehta
##
## Date Created: 11-03-2020
##
##
## ---------------------------
##
## Notes:
## Run the script with python3 unless some parts of the code will not work!
## This script is run with (ensure python3 is sourced): python3 ATAC_Bioinf_pipeline_v2b_part3a.py blastnoutput.blast mtgenome.fasta e.g. python3 ATAC_Bioinf_pipeline_v2b_part3a.py M_zebra_v1.1.assembly_mitochondrialBLASTn.blast Mz_mtgenome.fasta
## a. takes raw blast output - ensure that output is *.blast
## b. if file is empty, then confirm and echo as no mtDNA scaffolds for spID
## c. if pident=100 and evalue=0 then confirm sseqid as mtDNA scaffold
## d. if
##   pident>93; pidTresh=93, and
##   evalue<1e-10; eTresh=1e-10
##   then retain those hits
## Output file will be *filtered.blast
##
##
# # BLAST output headers are:
# 1. 	 qseqid 	 query (e.g., gene) sequence id
#  2. 	 sseqid 	 subject (e.g., reference genome) sequence id
#  3. 	 pident 	 percentage of identical matches
#  4. 	 length 	 alignment length
#  5. 	 mismatch 	 number of mismatches
#  6. 	 gapopen 	 number of gap openings
#  7. 	 qstart 	 start of alignment in query
#  8. 	 qend 	 end of alignment in query
#  9. 	 sstart 	 start of alignment in subject
#  10. 	 send 	 end of alignment in subject
#  11. 	 evalue 	 expect value
#  12. 	 bitscore 	 bit score
## ---------------------------

################################################################################################################

import sys

if len(sys.argv) < 2:
    print('\nNot enough arguments entered.\nUsage: python ATAC_Bioinf_pipeline_v2b_part3a.py blastn_output.blast mtgenome.fasta\n')
    sys.exit()

print("\n>>Importing libraries...")

import csv
import time
import collections
import os

print("\n>>Checking if BLAST input file is empty...\n")

## define a function to exit script if the input file is empty and echo 'BLAST output file is empty'
def is_file_empty(file_path):
    """ Check if blastn input file is empty by confirming if its size is 0 bytes"""
    # Check if file exist and it is empty
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0

# check if BLAST input file exists and it is not empty
is_empty = is_file_empty(sys.argv[1])
if is_empty:
    print('No genome scaffolds matching mitochondrial genome: BLAST output file is empty')
    print('...EXITING')
else:
    print('\nSome genome scaffolds match the mitochondrial genome')
    print('----------\nWorking with BLAST results from:\t' + sys.argv[1])
    print('----------\nWorking with mitochondrial genome sequence from:\t' + sys.argv[2])
    print('\n\n>>Filtering low quality hits...')

    # with open('M_zebra_v1.1.assembly_mitochondrialBLASTn.txt', 'r') as f:
    with open(sys.argv[1], 'r') as f:
        blast_out = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

    filt_dict=[]
    pidTresh=93
    eTresh=1e-10
    for x in blast_out:
        if float(x[2]) >= pidTresh:
            if float(x[10]) <= float(eTresh): # work with as a float as there are string values e.g. 1e-10
                filt_dict.append(x)
            else:
                filt_dict.append('no hits')
            # print(filt_dict)

    print('\nBLAST filtering done\n----------')

    ### At this stage, only accept matches that are hitting more than 75% of the full sequence length of the query (mitochondrial genome)
    # 1. take the length of the mitochondrial genome and store onto BLAST result
    # 2. take the alignment length column (3) and calculate a percentage of alignment length/total mtgenome length
    # 3. add the percentage as last column and only output those with >75%


    # This stores the FASTA file as a dictionary where the ID is key and whole sequence (not each line) is one value
    name = None
    seqs = dict()
    with open(sys.argv[2], 'r') as ref_seqs:
        for line in ref_seqs:
            line = line.rstrip() #discard the newline at the end (if any)
            #distinguish header from sequence
            if line.startswith('>'): #or line[0]=='>'
                #it is the header
                name = line[1:] #discarding the initial >
                seqs[name] = ''
            else:
                #it is sequence
                seqs[name] = seqs[name] + line

    print('\n\n>>Mitochondrial genome stored as dictionary and filtering based BLAST based on percentage hits to mt genome...')

    out_file = open(sys.argv[1][:-6] + '.filtered.blast', 'w') # open a file for writing out
    # out_file = open('M_zebra_v1.1.assembly_mitochondrialBLASTn.txt'[:-3] + 'filtered.blast', 'w') # open a file for writing out

    for id,seq in seqs.items():
        # print("mitochondrial FASTA ID is:\t" + id)
        # print("mitochondrial seqence is:\n" + seq)
        for hit in filt_dict: # take each line from the filtered BLAST results above
            # print(hit)
            # print(hit[0])
            # print("BLAST output lines printed")
            if hit[0] == id:
                res_list    = [] # create an empty list
                ref_seq_len = len(seq) # calculate the length of the sequence
                perc = round(int(hit[3])/int(ref_seq_len)*100) # calculate the percentage of alignment to the sequence length
                # print(hit,ref_seq_len,perc) # this prints the BLAST output, length of mtDNA match and (rounded) percentage aligned
                if perc >= 75: # if the % aligned is >=75% then store in res_list
                    res_list.append(hit + [ref_seq_len] + [perc]) # this appends ref_seq_len and perc to the hits list (and not appends onto end)
                    # print(res_list)
                    # print('\t'.join(str(x) for x in res_list[0]))
                    out_file.writelines('\t'.join(str(x) for x in res_list[0])) # as there are integers in the list, convert them to strings with for loop and print tab-delimited
            if hit[0] != id:
                # print("no hits")
                out_file.write("no hits") # can add a prompt to shell script to ignore files with "no hit" string and move on
    out_file.close()

    print('\nPercentage match filtering DONE and Outfile CLOSED\n----------')

    ######################################################
