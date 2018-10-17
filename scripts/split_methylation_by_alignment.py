#!/usr/bin/env python
# Split supplementary alignments from the output of nanopolish call-methylation
# Usage: python split_methylation_by_alignment.py /path/to/alignment.bam
# [/path/to/methylation.tsv]
from __future__ import print_function
import csv
import sys
import pickle
import copy


def get_read_id(read, contig, position, suppdb):
    if read in suppdb:
        for i, alignment in enumerate(suppdb[read]):
            chr, start, end = alignment
            if contig == chr and position >= start and position <= end:
                return "{}_{}".format(read, i)
        raise Exception("Position {}:{} not internal to any alignment of read {} ({})".format(
            contig, position, read, suppdb[read]))
    else:
        return read

bam_fn = sys.argv[1]
if len(sys.argv) > 2:
    meth_fn = sys.argv[2]
    handle = open(meth_fn, 'r')
else:
    print("Reading from stdin.", file=sys.stderr)
    handle = sys.stdin

with open("{}.suppdb".format(bam_fn), 'rb') as db_handle:
    suppdb = pickle.load(db_handle)

reader = csv.DictReader(handle, delimiter="\t")
fieldnames = reader.fieldnames
if "strand" in reader.fieldnames:
    fieldnames = copy.copy(fieldnames)
    del fieldnames[fieldnames.index("strand")]
    del_strand = True
else:
    del_strand = False
writer = csv.DictWriter(sys.stdout, delimiter="\t",
                        fieldnames=fieldnames)
writer.writeheader()
for row in reader:
    row["read_name"] = get_read_id(
        row["read_name"], row["chromosome"], int(row["start"]), suppdb)
    if del_strand:
        del row["strand"]
    try:
        writer.writerow(row)
    except Exception:
        print(row, file=sys.stderr)
        raise
handle.close()
