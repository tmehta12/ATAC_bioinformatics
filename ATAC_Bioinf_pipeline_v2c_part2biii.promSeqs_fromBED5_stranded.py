#!/usr/bin/python
import sys

if len(sys.argv) < 3:
    print('\nNot enough arguments entered...\n')
    print('Usage: genome FASTA, 5kb promoter BED\n')
    sys.exit()

print("\nImporting...\n")
import os
import csv
# import re
import operator
import itertools
import Bio
import time

from operator import itemgetter
from itertools import groupby
from Bio import SeqIO
from tqdm import tqdm

print('Loading\t' + str(sys.argv[1]) + '...')
records = SeqIO.parse(str(sys.argv[1]), "fasta")
print('Creating sequence dictionary...')
seq_dict = {}
for rec in tqdm(records):
    if rec.id in seq_dict:
        print('\n\tBreak due to error in seq dict constuction\n\n\n')
        sys.exit()
    else:
        seq_dict[rec.id] = rec.seq
print('Done.\n')

print('Working with annotations from:\t' + str(sys.argv[2] + '\n'))

# print('Writing to output file:\t' + str(sys.argv[2][:-12]) + 'fasta')
# out_file = open(sys.argv[2][:-12] + 'fasta', 'w')
out_file = open(sys.argv[2].replace('.bed', '.fasta'), 'w')
print('Writing to output file:\t' + sys.argv[2].replace('.bed', '.fasta'))
# sys.exit()

# for rec in records:
#     print(rec.id)
#
# for rec in records:
#     print ('Looking on chr:\t' + rec.id)
#     time.sleep(.1)
with open(sys.argv[2]) as bed_file:
    for gene in csv.reader(bed_file, delimiter='\t'):
        # print(gene)
        # print (gene[0])
        print('\n-----------------------------------------------------------\n')
        if gene[1] != gene[2]:
            if gene[0] in seq_dict:
                print ('~~~ On chromosome\t' + gene[0])
                print ('~~~ Retreiving gene\t' + gene[3])
                print ('~~~ From strand\t\t' + gene[5] + '\n')
                    # time.sleep(.1)
                prom_start = gene[1]
                print ('~~~ Promoter start pos:\t' + prom_start)
                prom_stop  = str(int(gene[2]))
                print ('~~~ Promoter stop pos:\t'  + prom_stop + '\n')

                # print(len(rec.seq))
                # print ('Sequence:\t' + rec.seq[int(prom_start):int(prom_stop)])

                if int(prom_stop) <= len(seq_dict[gene[0]]):
                    # print(gene[3])
                    out_file.write('>' + gene[3] + '\n')
                    if gene[5] == '+':
                        out_file.write(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)] + '\n'))
                        continue
                    elif gene[5] == '-':
                        out_file.write(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)].reverse_complement() + '\n'))
                        # print(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)].reverse_complement() + '\n'))
                        # print(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)].reverse_complement()[::-1] + '\n'))
                        continue
                else:
                    print('\n\nGENE START ' + str(prom_start) + ' at end of scaffold (scaffold end pos: ' + str(len(seq_dict[gene[0]])) + ')\n')

                    if gene[5] == '-' and len(seq_dict[gene[0]]) - int(prom_start) > 0:
                        prom_stop = str(len(seq_dict[gene[0]]))
                        print ('~~~ Revised promoter stop pos:\t'  + prom_stop + '\n')
                        out_file.write('>' + gene[3] + '\n')
                        # print(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)].reverse_complement()))
                        out_file.write(str(seq_dict[gene[0]][int(prom_start):int(prom_stop)].reverse_complement() + '\n'))
                        continue
                    else:
                        print('\n\t!!! NO LENGTH PROMOTER - Skipping...\n\n')

                # except:
                #     print('\nCould not write sequence for ' + gene[4])
                #     sys.exit()
        else:
            print ('~~~ On chromosome\t' + gene[0])
            print ('~~~ Retreiving gene\t' + gene[3])
            print ('~~~ Strand\t\t' + gene[5] + '\n')
                # time.sleep(.1)
            prom_start = gene[1]
            print ('~~~ Promoter start pos:\t' + prom_start)
            prom_stop  = str(int(gene[2]))
            print ('~~~ Promoter stop pos:\t'  + prom_stop + '\n')
            print('\n\t!!! NO LENGTH PROMOTER - Skipping...\n\n')
            # time.sleep(5)
            continue
    print('\n-----------------------------------------------------------\n')


out_file.close()
print('\nOutput file closed.\n\n')
