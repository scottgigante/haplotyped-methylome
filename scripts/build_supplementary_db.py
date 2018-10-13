#!/usr/bn/env python
# Builds a database of supplementary alignments in a bam file relative to read names
# Usage: ./build_supplementary_db.py /path/to/alignment.bam --summary out.tsv
from __future__ import print_function
import pysam
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bam')
parser.add_argument('-s', '--summary', default=None,
                    help='Optional file for printing read summaries')
parser.add_argument('--summary-only', default=False,
                    action='store_true', help="Don't produce a .suppdb file")
args = parser.parse_args()

bam_fn = args.bam
reads = dict()
with pysam.AlignmentFile(bam_fn, 'r') as bam:
    if args.summary:
        import numpy as np
        quals = dict()
    for read in bam:
        if read.is_unmapped:
            continue
        try:
            overlap = False
            for i in range(len(reads[read.query_name])):
                chr, start, end = reads[read.query_name][i]
                if chr == read.reference_name and start <= read.reference_end and end >= read.reference_start:
                    reads[read.query_name][i] = (
                        chr, min(start, read.reference_start), max(end, read.reference_end))
                    overlap = True
            if not overlap:
                reads[read.query_name].append(
                    (read.reference_name, read.reference_start, read.reference_end))
        except KeyError:
            reads[read.query_name] = [
                (read.reference_name, read.reference_start, read.reference_end)]
        if args.summary:
            try:
                read_qual = np.median(read.query_qualities)
            except TypeError:
                read_qual = "NA"
            try:
                quals[read.query_name].append(read_qual)
            except KeyError:
                quals[read.query_name] = [read_qual]
non_supplementary = set()
for read, alignments in reads.items():
    if len(alignments) == 1:
        non_supplementary.add(read)
if args.summary:
    with open(args.summary, 'w') as summary:
        print("read_name\tchr\tstart\tend\tqual", file=summary)
        for read, alignments in reads.items():
            if read in non_supplementary:
                print("{read_name}\t{chr}\t{start}\t{end}\t{qual}".format(
                    read_name=read,
                    chr=alignments[0][0],
                    start=alignments[0][1],
                    end=alignments[0][2],
                    qual=quals[read][0]), file=summary)
            else:
                for i, alignment in enumerate(alignments):
                    print("{read_name}\t{chr}\t{start}\t{end}\t{qual}".format(
                        read_name=read + "_" + str(i),
                        chr=alignment[0],
                        start=alignment[1],
                        end=alignment[2],
                        qual=quals[read][i]), file=summary)
if not args.summary_only:
    for read in non_supplementary:
        del reads[read]
    with open("{}.suppdb".format(bam_fn), 'wb') as handle:
        pickle.dump(reads, handle, pickle.HIGHEST_PROTOCOL)
