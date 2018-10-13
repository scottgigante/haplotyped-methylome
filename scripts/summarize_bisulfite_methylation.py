#!/usr/bin/env python
# Summarises bisulfite data per-site
# Input: Bismarck txt or txt.gz files with one methylation call (single read single site) per row
# Output: tsv with one line per genomic location
# Usage: ./summarize_bisulfite_methylation.py output.tsv input1.txt.gz
# [input2.txt.gz ...]

import sys
import csv
import gzip
__true = "+"
__false = "-"


def add(summary, chr, pos, meth):
    try:
        summary[chr]
    except KeyError:
        summary[chr] = dict()
    try:
        summary[chr][pos]
    except KeyError:
        summary[chr][pos] = [0., 0.]
    if meth:
        summary[chr][pos][0] += 1
    summary[chr][pos][1] += 1
    return summary


def summarize(in_handle, output_fn, summary=None):
    if summary is None:
        summary = dict()
    reader = csv.reader(in_handle, delimiter="\t")
#    header = next(reader)
    for line in reader:
        meth = line[1] == __true
        chr = line[2]
        pos = int(line[3])
        summary = add(summary, chr, pos, meth)
    return summary


def write_out_summary(output_fn, summary):
    with open(output_fn, 'w') as out_handle:
        for chr, positions in summary.items():
            for pos, data in positions.items():
                meth, total = tuple(data)
                out_handle.write("{}\t{}\t{:.3f}\t{}\t{}\n".format(
                    chr, pos, meth / total, meth, total))


output_fn = sys.argv[1]
input_fns = sys.argv[2:]
if len(sys.argv) == 3:
    input_fns = [sys.argv[2]]
else:
    input_fns = sys.argv[2:]

summary = None
for input_fn in input_fns:
    if input_fn.endswith("gz"):
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(input_fn, 'r') as handle:
        summary = summarize(handle, output_fn, summary)

write_out_summary(output_fn, summary)
