#!/usr/bin/env python

# Get read length (as integer) from fastq for TSS enrichment calculation
# Usage: python3

import os
import gzip
import re
import subprocess

# define the function
def get_read_length(fastq):
    # code adapted from Daniel Kim's ATAQC module
    def getFileHandle(filename, mode="r"):
        if (re.search('.gz$', filename) or re.search('.gzip', filename)):
            if (mode == "r"):
                mode = "rb"
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1
    return int(max_length)

read_len = get_read_length(sys.argv[1])

out_file = open(sys.argv[1] + '_read_length.txt', 'w')
out_file.write(read_len)
out_file.close()
