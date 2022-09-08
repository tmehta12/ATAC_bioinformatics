#!/usr/bin/python

#=====================================================================================================================================================

# 1. Adjusted peak coords according to alignment gaps using script 'peak_promaln_getcoords.py' and ran in 'ATAC_Bioinf_ManuscriptAnalysis.sh'
# 2. Each orthogroups new bed files were joined by tissue in 'ATAC_Bioinf_ManuscriptAnalysis.sh')
# 3. In 'ATAC_Bioinf_ManuscriptAnalysis.sh', bedtools intersect was ran to determine peaks that overlap between species - reporting both ways of overlap so Ab>Mz and Mz>Ab despite redundancy
# 4. In 'ATAC_Bioinf_ManuscriptAnalysis.sh', then checked that both species peak summits are within the overlapping region

# Summary: This script will then calculate % pairwise identity over the overlapping peaks between species to study genetic diversity

# In steps:
# - read_msa : alignment fasta to ordered dict (reads a multiple a sequence alignment file in FASTA/Pearson format)
# - read_bed : bed file as list of lists
# - compare : loop through bed file and for each pairwise overlap, compute percentage sequence similarity
    # a. read in the fasta as an ordered dict
    # b. read in the bed file with overlapping peaks and store each column as list of lists
    # c. store the ref and alt sequences for each bed line
    # d. using the range of peak_overlap_start and peak_overlap_end, calculate and store identity, similarity and their percentages (excluding gaps)
    # e. output as a list of lists
# output the original input bed file with additional cols (alignment_length, identical residues and percentage identity)

# output file format:
# 1.refensID_bed
# 2.refpeakaln_start
# 3.refpeakaln_end
# 4.refpeakID
# 5.refscore
# 6.refstrand
# 7.refMACS2_signalValue
# 8.refMACS2_pval
# 9.refMACS2_qval
# 10.refpeakaln_summit
# 11.refIDR
# 12.reforthogroup
# 13.refens_gs1
# 14.refens_gs2
# 15.refprom_chr
# 16.refprom_start
# 17.refprom_end
# 18.refpeaklength
# 19.reftss_coord
# 20.refpeak_genchr
# 21.refpeak_genstart
# 22.refpeak_genend
# 23.refpeak_start
# 24.refpeak_end
# 25.refpeak_summit
# 26.refspID
# 27.altensID_bed
# 28.altpeakaln_start
# 29.altpeakaln_end
# 30.altpeakID
# 31.altscore
# 32.altstrand
# 33.altMACS2_signalValue
# 34.altMACS2_pval
# 35.altMACS2_qval
# 36.altpeakaln_summit
# 37.altIDR
# 38.altorthogroup
# 39.altens_gs1
# 40.altens_gs2
# 41.altprom_chr
# 42.altprom_start
# 43.altprom_end
# 44.altpeaklength
# 45.alttss_coord
# 46.altpeak_genchr
# 47.altpeak_genstart
# 48.altpeak_genend
# 49.altpeak_start
# 50.altpeak_end
# 51.altpeak_summit
# 52.altspID
# 53.peakoverlap
# 54.peak_overlap_start
# 55.peak_overlap_end
# 56.alignment_length (excludes gaps)
# 57.identical residues
# 58.percentage identity

#=====================================================================================================================================================

print('\n===================================================================================\n')
print('peak_promaln_overlaps_gd.py')
import sys

if len(sys.argv) < 3:
    print('\nNot enough arguments entered.\nUsage: python3 peak_promaln_overlaps_gd.py alignmentFASTA.fa peakalnoverlapBEDfile outfile\n')
    print('\n===================================================================================\n')
    sys.exit()


print("\nImporting...\n\n")
import os
import csv
import collections
import itertools
import re

from collections import OrderedDict

# 1. input files
fasta_file = sys.argv[1]
# fasta_file = 'OMA00020_aln.fa'

bed_file = sys.argv[2]
# bed_file = 'OMA00020.proalnmap.peakovrlap.SP.bed'

# 2. function to read in the FASTA alignment file (in FASTA/Pearson format) as an ordered dictionary
def read_msa(file_name):
    seq_dict = OrderedDict()
    try:
        with open(file_name, "r") as seqfile:
            for line in seqfile:
                line = line.strip()
                if line[0] == ">":
                    seq_dict[line[1:]] = ""
                else:
                    seq_dict[next(reversed(seq_dict))] = seq_dict[next(reversed(seq_dict))] + line
    except:
        print(fasta_file + " could not be found, exiting!")
    return seq_dict

# 3. function to open the peak overlap bed file and for each line, strip by newline, split by tab and assign to bed_out list
def read_bed(bed_file_name):
    try:
        with open(bed_file_name, "r") as bed_file_in:
            bed_list = [line.strip('\n').split('\t') for line in bed_file_in]
    except:
        print(bed_file + " could not be found, exiting!")
    return bed_list


# 3. get the percentage genetic diversity over the overlapping region
# a. read in the fasta as an ordered dict
# b. read in the bed file with overlapping peaks and store each column as list of lists
# c. store the ref and alt sequences for each bed line
# d. using the range of peak_overlap_start and peak_overlap_end, calculate and store identity, similarity and their percentages (excluding gaps)
# e. output as a list of lists

seq_dict = read_msa(fasta_file)
bed_list = read_bed(bed_file)
bed_out = [] # bed out as list

def compare():
    for row in bed_list: # map each col in the BED to a variable
        refensID_bed = row[0]
        refpeakaln_start = row[1]
        refpeakaln_end = row[2]
        refpeakID = row[3]
        refscore = row[4]
        refstrand = row[5]
        refMACS2_signalValue = row[6]
        refMACS2_pval = row[7]
        refMACS2_qval = row[8]
        refpeakaln_summit = row[9]
        refIDR = row[10]
        reforthogroup = row[11]
        refens_gs1 = row[12]
        refens_gs2 = row[13]
        refprom_chr = row[14]
        refprom_start = row[15]
        refprom_end = row[16]
        refpeaklength = row[17]
        reftss_coord = row[18]
        refpeak_genchr = row[19]
        refpeak_genstart = row[20]
        refpeak_genend = row[21]
        refpeak_start = row[22]
        refpeak_end = row[23]
        refpeak_summit = row[24]
        refspID = re.sub(r'[0-9]+', '', refpeakID.split('_')[0]) # create a species ID variable from the peakID
        altensID_bed = row[25]
        altpeakaln_start = row[26]
        altpeakaln_end = row[27]
        altpeakID = row[28]
        altscore = row[29]
        altstrand = row[30]
        altMACS2_signalValue = row[31]
        altMACS2_pval = row[32]
        altMACS2_qval = row[33]
        altpeakaln_summit = row[34]
        altIDR = row[35]
        altorthogroup = row[36]
        altens_gs1 = row[37]
        altens_gs2 = row[38]
        altprom_chr = row[39]
        altprom_start = row[40]
        altprom_end = row[41]
        altpeaklength = row[42]
        alttss_coord = row[43]
        altpeak_genchr = row[44]
        altpeak_genstart = row[45]
        altpeak_genend = row[46]
        altpeak_start = row[47]
        altpeak_end = row[48]
        altpeak_summit = row[49]
        altspID = re.sub(r'[0-9]+', '', altpeakID.split('_')[0]) # create a species ID variable from the peakID
        peakoverlap = row[50]
        peak_overlap_start = row[51]
        peak_overlap_end = row[52]
        iden = 0
        sim = 0
        seq1 = seq_dict[refensID_bed]
        seq2 = seq_dict[altensID_bed]
        align_len = 0
        # min_len = min(len(seq1), len(seq2))
        for i in range(int(peak_overlap_start),int(peak_overlap_end)):
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
        align_length = str(align_len)
        identity = str(iden)
        perc_iden = str(round(iden / align_len * 100,2))
        bed_out.extend([[refensID_bed,refpeakaln_start,refpeakaln_end,refpeakID,refscore,refstrand,refMACS2_signalValue,refMACS2_pval,refMACS2_qval,refpeakaln_summit,refIDR,reforthogroup,refens_gs1,refens_gs2,refprom_chr,refprom_start,refprom_end,refpeaklength,reftss_coord,refpeak_genchr,refpeak_genstart,refpeak_genend,refpeak_start,refpeak_end,refpeak_summit,refspID,altensID_bed,altpeakaln_start,altpeakaln_end,altpeakID,altscore,altstrand,altMACS2_signalValue,altMACS2_pval,altMACS2_qval,altpeakaln_summit,altIDR,altorthogroup,altens_gs1,altens_gs2,altprom_chr,altprom_start,altprom_end,altpeaklength,alttss_coord,altpeak_genchr,altpeak_genstart,altpeak_genend,altpeak_start,altpeak_end,altpeak_summit,altspID,peakoverlap,peak_overlap_start,peak_overlap_end,align_length,identity,perc_iden]])

compare() # run the function to generate the output

with open(sys.argv[3], 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in bed_out) # join each element as tab-separated and trailing newlines for each line in the list
file.close()

print('\nOut bed file written with file name:\t' + str(sys.argv[3]))
