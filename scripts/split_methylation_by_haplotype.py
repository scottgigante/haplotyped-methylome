#!/usr/bin/env python
# Usage: python compare_haplotype_methylation.py methylation.split_supplementary.sorted.per_site.tsv haplotypes.tsv
# Compares methylation states between groups of haplotyped methylation calls along the genome using a rolling t-test
# Haplotypes file should have one line per read (where supplementary
# alignments are counted separately). Methylation file should have one
# line per read per site, sorted by site.

from __future__ import print_function
import sys
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--methylation', required=True,
                    help="Methylation data, split by supplementaries, sorted by genomic location")
parser.add_argument('-p', '--phase', required=True,
                    help="Phased reads (haplotyping data)")
args = parser.parse_args()

meth_fn = args.methylation
haplotype_fn = args.phase

haplotypes = dict()
genotypes = set(["fail"])
with open(haplotype_fn, 'r') as haplotype_handle:
    reader = csv.DictReader(haplotype_handle, delimiter="\t")
    for row in reader:
        haplotypes[row["read"]] = row["genotype"]
        genotypes.add(row["genotype"])


def get_genotype(read_name):
    try:
        return haplotypes[read_name]
    except KeyError:
        return "fail"

with open(meth_fn, 'r') as meth_handle:
    reader = csv.DictReader(meth_handle, delimiter="\t")
    out_handles = list()
    writers = dict()
    for genotype in genotypes:
        handle = open("{}.{}.tsv".format(meth_fn, genotype), 'w')
        out_handles.append(handle)
        writers[genotype] = csv.DictWriter(
            handle, fieldnames=reader.fieldnames, delimiter="\t")
        writers[genotype].writeheader()
    for row in reader:
        writers[get_genotype(row["read_name"])].writerow(row)
    for handle in out_handles:
        handle.close()
