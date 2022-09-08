#!/usr/bin/python
#=====================================================================================================================================================

# for every pairwise instance that DOES NOT exist e.g. Mz-On, Mz-Nb and Mz-Pn present but Mz-Ab absent, gets marked as a loss in Ab
# read each line comparing the pairwise presence (Mz, Pnm, Ab, Nb only) for a ref peak, and then output:
    # 1. orthogroup
    # 2. ref peak
    # 3. comparison species peak
    # 4. genetic diversity
    # 5. species peak presence
    # 6. losses in other species
# Input:
    # OMA19993:On-On3_B_ATAC_peak_111231b 	Mz-Mz2_B_ATAC_peak_104899b:63.12	Nb-Nb4_B_ATAC_peak_136281c:64.57	Nb-Nb5_B_ATAC_peak_146910b:64.57
# Output:
    # OMA19993  On3_B_ATAC_peak_111231b  Mz2_B_ATAC_peak_104899b  63.12 PRESENCE:On-Mz  LOSS:On-Pnm;On-Ab
    # OMA19993  On3_B_ATAC_peak_111231b  Nb4_B_ATAC_peak_136281c  64.57 PRESENCE:On-Mz  LOSS:On-Pnm;On-Ab
    # OMA19993  On3_B_ATAC_peak_111231b  Nb5_B_ATAC_peak_146910b  64.57 PRESENCE:On-Mz  LOSS:On-Pnm;On-Ab

#=====================================================================================================================================================

print('\n===================================================================================\n')
print('peak_promaln_gdloss.py')
import sys

if len(sys.argv) < 2:
    print('\nNot enough arguments entered.\nUsage: python3 peak_promaln_gdloss.py peak_file outfile\n')
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
peak_file = sys.argv[1]
# peak_file = 'B.proalnmap.peakovrlap.SP.gd.collated.1to1.loss_2.txt'

# 2. function to read the tab-delimited peak file, line by line, and create a dictionary (key is first col) and other cols are values
def read_peaks(file_name):
    try:
        with open(file_name, 'r') as peakfile:
             rows = (line.strip('\n').split('\t') for line in peakfile)
             dict = {row[0]:row[1:] for row in rows}
    except:
        print(peak_file + " could not be found, exiting!")
    return dict

# 3. Determine the species specific peak losses
# a. read in the peak file as dictionary
# b. read dict line by line
# c. then scan and store the losses by comparing against the presence in the altpeaks values
# d. extract the ref species, peaks, all alt species, peaks and genetic diversity

peak_dict = read_peaks(peak_file)
loss_list = [] # create an empty list to store the losses

for orthrefpeak, altpeaks in peak_dict.items():
    # print(orthrefpeak, ' : ', altpeaks) # access the keys and values
    sp1 = 'Mz'
    sp2 = 'Pnm'
    sp3 = 'Ab'
    sp4 = 'Nb'
    sp5 = 'On'
    orthogroup = orthrefpeak.split(':')[0]
    refsppeak = orthrefpeak.split(':')[1]
    refsp = refsppeak.split('-')[0]
    refpeak = refsppeak.split('-')[1]
    # print(orthogroup,refsp,refpeak,altpeaks)
    for x in altpeaks:
        altsp = x.split('-')[0]
        altpeak = x.split('-')[1].split(':')[0]
        geneticdiversity = x.split(':')[1]
        # print(orthogroup,refpeak,altpeak,geneticdiversity,'PRESENCE:'+refsp+'-'+altsp)
        # print(orthogroup,refpeak,altpeak,geneticdiversity,'PRESENCE:'+refsp+'-'+altsp,altpeaks)
        if sp1 in orthrefpeak and sp2 not in '\t'.join(altpeaks):
            # print(orthogroup,refpeak,altpeak,geneticdiversity,'PRESENCE:'+refsp+'-'+altsp,'LOSS:'+sp1+'-'+sp2)
            comp1 = 'LOSS:'+sp1+'-'+sp2+';'
        else:
            # print(orthogroup,refpeak,altpeak,geneticdiversity,'PRESENCE:'+refsp+'-'+altsp)
            comp1 = 'LOSS:NA-NA;'
        if sp1 in orthrefpeak and sp3 not in '\t'.join(altpeaks):
            comp2 = 'LOSS:'+sp1+'-'+sp3+';'
        else:
            comp2 = 'LOSS:NA-NA;'
        if sp1 in orthrefpeak and sp4 not in '\t'.join(altpeaks):
            comp3 = 'LOSS:'+sp1+'-'+sp4+';'
        else:
            comp3 = 'LOSS:NA-NA;'
        if sp2 in orthrefpeak and sp1 not in '\t'.join(altpeaks):
            comp4 = 'LOSS:'+sp2+'-'+sp1+';'
        else:
            comp4 = 'LOSS:NA-NA;'
        if sp2 in orthrefpeak and sp3 not in '\t'.join(altpeaks):
            comp5 = 'LOSS:'+sp2+'-'+sp3+';'
        else:
            comp5 = 'LOSS:NA-NA;'
        if sp2 in orthrefpeak and sp4 not in '\t'.join(altpeaks):
            comp6 = 'LOSS:'+sp2+'-'+sp4+';'
        else:
            comp6 = 'LOSS:NA-NA;'
        if sp3 in orthrefpeak and sp1 not in '\t'.join(altpeaks):
            comp7 = 'LOSS:'+sp3+'-'+sp1+';'
        else:
            comp7 = 'LOSS:NA-NA;'
        if sp3 in orthrefpeak and sp2 not in '\t'.join(altpeaks):
            comp8 = 'LOSS:'+sp3+'-'+sp2+';'
        else:
            comp8 = 'LOSS:NA-NA;'
        if sp3 in orthrefpeak and sp4 not in '\t'.join(altpeaks):
            comp9 = 'LOSS:'+sp3+'-'+sp4+';'
        else:
            comp9 = 'LOSS:NA-NA;'
        if sp4 in orthrefpeak and sp1 not in '\t'.join(altpeaks):
            comp10 = 'LOSS:'+sp4+'-'+sp1+';'
        else:
            comp10 = 'LOSS:NA-NA;'
        if sp4 in orthrefpeak and sp2 not in '\t'.join(altpeaks):
            comp11 = 'LOSS:'+sp4+'-'+sp2+';'
        else:
            comp11 = 'LOSS:NA-NA;'
        if sp4 in orthrefpeak and sp3 not in '\t'.join(altpeaks):
            comp12 = 'LOSS:'+sp4+'-'+sp3+';'
        else:
            comp12 = 'LOSS:NA-NA;'
        if sp5 in orthrefpeak and sp1 not in '\t'.join(altpeaks):
            comp13 = 'LOSS:'+sp5+'-'+sp1+';'
        else:
            comp13 = 'LOSS:NA-NA;'
        if sp5 in orthrefpeak and sp2 not in '\t'.join(altpeaks):
            comp14 = 'LOSS:'+sp5+'-'+sp2+';'
        else:
            comp14 = 'LOSS:NA-NA;'
        if sp5 in orthrefpeak and sp3 not in '\t'.join(altpeaks):
            comp15 = 'LOSS:'+sp5+'-'+sp3+';'
        else:
            comp15 = 'LOSS:NA-NA;'
        if sp5 in orthrefpeak and sp4 not in '\t'.join(altpeaks):
            comp16 = 'LOSS:'+sp5+'-'+sp4
        else:
            comp16 = 'LOSS:NA-NA'
        loss_list.extend([[orthogroup,refpeak,altpeak,geneticdiversity,'PRESENCE:'+refsp+'-'+altsp,comp1+comp2+comp3+comp4+comp5+comp6+comp7+comp8+comp9+comp10+comp11+comp12+comp13+comp14+comp15+comp16]])

loss_list = [[item.replace('LOSS:NA-NA;', '') for item in lst] for lst in loss_list]
loss_list = [[item.replace(';LOSS:NA-NA', '') for item in lst] for lst in loss_list]

# for line in loss_list:
#     print(line)

with open(sys.argv[2], 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in loss_list) # join each element as tab-separated and trailing newlines for each line in the list
file.close()

print('\nOut bed file written with file name:\t' + str(sys.argv[2]))
