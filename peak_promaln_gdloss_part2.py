#!/usr/bin/python

#=====================================================================================================================================================

# 1. Adjusted peak coords according to alignment gaps using script 'peak_promaln_getcoords.py' and ran in 'ATAC_Bioinf_ManuscriptAnalysis.sh'
# 2. Each orthogroups new bed files were joined by tissue in 'ATAC_Bioinf_ManuscriptAnalysis.sh')
# 3. In 'ATAC_Bioinf_ManuscriptAnalysis.sh', bedtools intersect was ran to determine peaks that overlap between species - reporting both ways of overlap so Ab>Mz and Mz>Ab despite redundancy
# 4. In 'ATAC_Bioinf_ManuscriptAnalysis.sh', then checked that both species peak summits are within the overlapping region
# 5. Peak loss determined using python script: peak_promaln_gdloss.py
# 6. In 'ATAC_Bioinf_ManuscriptAnalysis.sh', the files were then prepped for this script

# Summary: This script will then calculate % pairwise identity over the absent peak regions (so loss) to study genetic diversity

# In steps:
# - read_peaks : peak file as list of lists
# - read_orth : orthologs file as dictionary
# - read_msa : alignment fasta to ordered dict (reads a multiple a sequence alignment file in FASTA/Pearson format)
# - compare : loop through peak file and for each pairwise species loss, compute percentage sequence similarity between the ref peak region and absent species region
    # a. read in the overlapping and loss/gain peaks file and store each column as list of lists
    # b. read in the orthologs file as a dictionary
    # c. from the peaks list orthogroup of each line, read in the corresponding alignment fasta file as an ordered dict
    # d. from the peaks list, determine if the line is a pairwise species loss and if so, then use the refpeakaln_start and refpeakaln_end, to calculate and store sequence similarity and their percentages (excluding gaps) between the ref peak region and the corresponding non peak region in the comparison species (with a peak loss)
    # e. output as a list of lists with all other lines (non loss or gain as NA for perc_iden)

# output the original input bed file with additional cols (percentage identity)

# output file format:
# 1.orthogroup
# 2.refpeakID
# 3.altpeakID
# 4.ref_alt_peak_geneticdiversity
# 5.ref_alt_peak_presencesp
# 6.ref_alt_peak_gainlosssp
# 7.refpeak_alt_loss_percentage_identity
# 8.reforthogroup
# 9.refpeakaln_start
# 10.refpeakaln_end
# 11.refpeakIDrep
# 12.refscore
# 13.refstrand
# 14.refMACS2_signalValue
# 15.refMACS2_pval
# 16.refMACS2_qval
# 17.refpeakaln_summit
# 18.refIDR
# 19.refensID_bed
# 20.refens_gs1
# 21.refens_gs2
# 22.refprom_chr
# 23.refprom_start
# 24.refprom_end
# 25.refpeaklength
# 26.reftss_coord
# 27.refpeak_genchr
# 28.refpeak_genstart
# 29.refpeak_genend
# 30.refpeak_start
# 31.refpeak_end
# 32.refpeak_summit

#=====================================================================================================================================================

print('\n===================================================================================\n')
print('peak_promaln_gdloss_part2.py')
import sys

if len(sys.argv) < 4:
    print('\nNot enough arguments entered.\nUsage: python3 peak_promaln_gdloss_part2.py [peak_aln_overlaploss_file] [path_to_alignment_fastas_with_trailingslash] [1to1orthfile] [outfile]\n')
    print('\n===================================================================================\n')
    sys.exit()


print("\nImporting...\n\n")
import os
import csv
import collections
import itertools
import re

from collections import OrderedDict
print("\nImported.\n\n")

print("\nLoading files and paths...\n\n")
# 1. input files
peak_file = sys.argv[1]
# peak_file = 'B.zeropeakoverlap_GAIN.peakovrlap_LOSS.1to1.txt2.TEST'

fasta_path = sys.argv[2]
# fasta_path = '/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/4.Peak_Orth/peakcomp_gdiversity/orthoGroup_fasta/alignments/'
# fasta_path = '/Users/mehtat/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/'

orth_file = sys.argv[3]
# orth_file = 'OrthologousGroups_ENS_GeneSymbols_5cich1to1.collated'
print("\nFiles and paths loaded.\n\n")


# 2. function to open the peak overlap/loss file and for each line, strip by newline, split by tab and assign to list
def read_peaks(peak_file_name):
    try:
        with open(peak_file_name, "r") as peak_file_in:
            peak_list = [line.strip('\n').split('\t') for line in peak_file_in]
    except:
        print(peak_file + " could not be found, exiting!")
    return peak_list

# 3. function to read in the FASTA alignment file (in FASTA/Pearson format) as an ordered dictionary
def read_msa(fasta_file_name):
    seq_dict = OrderedDict()
    try:
        with open(fasta_file_name, "r") as seqfile:
            for line in seqfile:
                line = line.strip()
                if line[0] == ">":
                    seq_dict[line[1:]] = ""
                else:
                    seq_dict[next(reversed(seq_dict))] = seq_dict[next(reversed(seq_dict))] + line
    except:
        print(orthogroup + "_aln.fa could not be found, exiting!")
    return seq_dict

# 4. function to read the tab-delimited 1to1 orthologs file, line by line, and create a dictionary (key is first col) and other cols are values
def read_orth(orth_file_name):
    try:
        with open(orth_file_name, 'r') as orthfile:
             rows = (line.strip('\n').split('\t') for line in orthfile)
             dict = {row[0]:row[1:] for row in rows}
    except:
        print(orth_file + " could not be found, exiting!")
    return dict

print("\nFunctions to read files created.\n\n")

# 5. get the percentage genetic diversity over the overlapping region
# a. read in the overlapping and loss/gain peaks file and store each column as list of lists
# b. read in the orthologs file as a dictionary
# c. from the peaks list orthogroup of each line, read in the corresponding alignment fasta file as an ordered dict
# d. from the peaks list, determine if the line is a pairwise species loss and if so, then use the refpeakaln_start and refpeakaln_end, to calculate and store sequence similarity and their percentages (excluding gaps) between the ref peak region and the corresponding non peak region in the comparison species (with a peak loss)
# e. output as a list of lists with all other lines (non loss or gain as NA for perc_iden)

peak_list = read_peaks(peak_file)
# for row in peak_list:
#     print(row)
orth_dict = read_orth(orth_file)
peak_out = [] # peak out as list

print("\nFiles stored as lists and dictionaries.\n\n")

def compare():
    for row in peak_list: # map each col in the peak file to a variable
        orthogroup = row[0]
        refpeakID = row[1]
        altpeakID = row[2]
        ref_alt_peak_geneticdiversity = row[3]
        ref_alt_peak_presencesp = row[4]
        ref_alt_peak_gainlosssp = row[5]
        reforthogroup = row[6]
        refpeakaln_start = row[7]
        refpeakaln_end = row[8]
        refpeakIDrep = row[9]
        refscore = row[10]
        refstrand = row[11]
        refMACS2_signalValue = row[12]
        refMACS2_pval = row[13]
        refMACS2_qval = row[14]
        refpeakaln_summit = row[15]
        refIDR = row[16]
        refensID_bed = row[17]
        refens_gs1 = row[18]
        refens_gs2 = row[19]
        refprom_chr = row[20]
        refprom_start = row[21]
        refprom_end = row[22]
        refpeaklength = row[23]
        reftss_coord = row[24]
        refpeak_genchr = row[25]
        refpeak_genstart = row[26]
        refpeak_genend = row[27]
        refpeak_start = row[28]
        refpeak_end = row[29]
        refpeak_summit = row[30]
        seq_dict = read_msa(fasta_path+orthogroup+"_aln.fa")
        gainorloss = ref_alt_peak_gainlosssp.split(':')[0]
        if gainorloss == "LOSS":
            gainorlossna = ref_alt_peak_gainlosssp.split(':')[1]
            # print(gainorloss,gainorlossna)
            if gainorlossna != "NA-NA":
                # print(gainorloss,gainorlossna)
                altspID = ref_alt_peak_gainlosssp.split('-')[1]
                altensID_bed = orth_dict[orthogroup+"_"+altspID][0]
                iden = 0
                sim = 0
                seq1 = seq_dict[refensID_bed]
                seq2 = seq_dict[altensID_bed]
                align_len = 0
                # print(gainorloss,gainorlossna, altspID, seq1, seq2)
                # min_len = min(len(seq1), len(seq2))
                for i in range(int(refpeakaln_start),int(refpeakaln_end)):
                    if seq1[i] == "-" and seq2[i] == "-":
                        pass
                    else:
                        align_len += 1
                        if seq1[i] == seq2[i]:
                            iden += 1
                            sim += 1
                # print(refensID_bed + " vs " + altensID_bed)
                # print("Alignment length: ", align_len)
                # print("Identical residues: ", iden)
                # print("Percent identity: ", round(iden / align_len * 100,2))
                # print("Similar residues: ", sim-iden)
                # print("Percent similarity: " , round(sim / align_len * 100,2))
                # print(" ")
                # align_length = str(align_len)
                identity = str(iden)
                perc_iden = str(round(iden / align_len * 100,2))
                peak_out.extend([[orthogroup,refpeakID,altpeakID,ref_alt_peak_geneticdiversity,ref_alt_peak_presencesp,ref_alt_peak_gainlosssp,perc_iden,reforthogroup,refpeakaln_start,refpeakaln_end,refpeakIDrep,refscore,refstrand,refMACS2_signalValue,refMACS2_pval,refMACS2_qval,refpeakaln_summit,refIDR,refensID_bed,refens_gs1,refens_gs2,refprom_chr,refprom_start,refprom_end,refpeaklength,reftss_coord,refpeak_genchr,refpeak_genstart,refpeak_genend,refpeak_start,refpeak_end,refpeak_summit]])
            else:
                peak_out.extend([[orthogroup,refpeakID,altpeakID,ref_alt_peak_geneticdiversity,ref_alt_peak_presencesp,ref_alt_peak_gainlosssp,"NA",reforthogroup,refpeakaln_start,refpeakaln_end,refpeakIDrep,refscore,refstrand,refMACS2_signalValue,refMACS2_pval,refMACS2_qval,refpeakaln_summit,refIDR,refensID_bed,refens_gs1,refens_gs2,refprom_chr,refprom_start,refprom_end,refpeaklength,reftss_coord,refpeak_genchr,refpeak_genstart,refpeak_genend,refpeak_start,refpeak_end,refpeak_summit]]) # print the lines that have LOSS:NA-NA
        else:
            peak_out.extend([[orthogroup,refpeakID,altpeakID,ref_alt_peak_geneticdiversity,ref_alt_peak_presencesp,ref_alt_peak_gainlosssp,"NA",reforthogroup,refpeakaln_start,refpeakaln_end,refpeakIDrep,refscore,refstrand,refMACS2_signalValue,refMACS2_pval,refMACS2_qval,refpeakaln_summit,refIDR,refensID_bed,refens_gs1,refens_gs2,refprom_chr,refprom_start,refprom_end,refpeaklength,reftss_coord,refpeak_genchr,refpeak_genstart,refpeak_genend,refpeak_start,refpeak_end,refpeak_summit]]) # print the GAIN lines

print("\nFunction to compare created.\n\n")

print("\nRunning comparison....\n\n")

compare() # run the function to generate the output

# for row in peak_out:
#     print(row)

print("\nComparison completed and now writing file....\n\n")

with open(sys.argv[4], 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in peak_out) # join each element as tab-separated and trailing newlines for each line in the list
file.close()

# with open('B.test.out', 'w') as file: # open a file for writing
#     file.writelines('\t'.join(i) + '\n' for i in peak_out) # join each element as tab-separated and trailing newlines for each line in the list
# file.close()

print('\nOut bed file written with file name:\t' + str(sys.argv[4]))
