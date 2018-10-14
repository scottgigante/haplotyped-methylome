#!/usr/bin/env python
# Usage: python compare_haplotype_methylation.py methylation.split_supplementary.sorted.per_site.tsv haplotypes.tsv
# Compares methylation states between groups of haplotyped methylation calls along the genome using a rolling t-test
# Haplotypes file should have one line per read (where supplementary
# alignments are counted separately). Methylation file should have one
# line per read per site, sorted by site.

from __future__ import print_function
import csv
import numpy as np
import math
import argparse
import re


def parse_region(region):
    """
    Parse a region specification.
    :param region: String region specification
    :raises argparse.ArgumentTypeError: raises an error if format not recognised
    :returns contig: String contig / chromosome name
    :returns start: integer start position (0-based)
    :returns end: integer end position (1-based)
    >>> parseRegion("chr1:1000-2000")
    ("chr1", 1000, 2000)
    """
    r = ''.join(region.split())  # remove whitespace
    r = re.split(':|-', r)
    start = 0
    end = None
    if len(r) < 1:
        raise argparse.ArgumentTypeError(
            "Region {} must specify a reference name".format(region))
    elif len(r) > 3:
        raise argparse.ArgumentTypeError(
            "Region {} format not recognized".format(region))
    else:
        contig = r[0]
        try:
            start = int(re.sub(",|\.", "", r[1]))
        except IndexError:
            pass
        try:
            end = int(re.sub(",|\.", "", r[2]))
        except IndexError:
            pass
        finally:
            return contig, start, end


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--methylation', required=True,
                    help="Methylation data, split by supplementaries, sorted by genomic location")
parser.add_argument('-p', '--phase', required=True,
                    help="Phased reads (haplotyping data)")
parser.add_argument('-b', '--binwidth', type=int, default=11,
                    help="Genomic distance over which to aggregate tests")
parser.add_argument('-o', '--overlap', type=int, default=5,
                    help="Genomic distance for bin overlap. Limited to binwidth/2 for streaming reasons, but this could be reengineered")
parser.add_argument('-r', '--region', type=parse_region, nargs="+", metavar="CHR:START-END",
                    help="Name and range of contigs to test (Default: entire genome)")
args = parser.parse_args()

__binwidth = args.binwidth
__overlap = min(__binwidth // 2, args.overlap)  # max = __binwidth/2

meth_fn = args.methylation
haplotype_fn = args.phase

regions = dict()
# format: { chr : [start, end] }
if args.region:
    for region in args.region:
        try:
            regions[region[0]]
        except KeyError:
            regions[region[0]] = []
        regions[region[0]].append([region[1], region[2]])

haplotypes = dict()
with open(haplotype_fn, 'r') as haplotype_handle:
    reader = csv.DictReader(haplotype_handle, delimiter="\t")
    for row in reader:
        haplotypes[row["read"]] = row["genotype"]


def get_genotype(read_name):
    try:
        return haplotypes[read_name]
    except KeyError:
        return "fail"


print("\t".join(["chr", "start", "end", "ref_mean", "ref_sd", "ref_site_count", "ref_read_count",
                 "alt_mean", "alt_sd", "alt_site_count", "alt_read_count"]))
chr = None
start = None
ref_reads = None
alt_reads = None
ref_meth = None
alt_meth = None
next_start = None
next_ref_reads = None
next_alt_reads = None
next_ref_meth = None
next_alt_meth = None
end = 0

with open(meth_fn, 'r') as meth_handle:
    reader = csv.DictReader(meth_handle, delimiter="\t")
    for row in reader:
        pos = (float(row["start"]) + float(row["end"])) / 2
        if regions:
            hit = False
            try:
                for region in regions[row['chromosome']]:
                    if pos >= region[0] and pos < region[1]:
                        hit = True
            except KeyError:
                pass
            if not hit:
                continue
        if pos >= end or row["chromosome"] != chr:
            if chr is not None:
                # print result
                if not (len(ref_meth) == 0 and len(alt_meth) == 0):
                    ref_mean = np.mean(ref_meth) if len(ref_meth) > 0 else "NA"
                    ref_sd = np.std(ref_meth) if len(ref_meth) > 0 else "NA"
                    alt_mean = np.mean(alt_meth) if len(alt_meth) > 0 else "NA"
                    alt_sd = np.std(alt_meth) if len(alt_meth) > 0 else "NA"
                    print("\t".join([chr, str(start), str(end), str(ref_mean), str(ref_sd), str(len(ref_meth)), str(len(
                        ref_reads)), str(alt_mean), str(alt_sd), str(len(alt_meth)), str(len(alt_reads))]))
            if chr is not None and chr == row["chromosome"] and pos < next_start:
                start = next_start
                ref_reads = next_ref_reads
                alt_reads = next_alt_reads
                ref_meth = next_ref_meth
                alt_meth = next_alt_meth
            else:
                start = int(row["start"])
                ref_reads = set()
                alt_reads = set()
                ref_meth = list()
                alt_meth = list()
                chr = row["chromosome"]
            next_start = start + __overlap
            end = start + __binwidth
            next_ref_reads = set()
            next_alt_reads = set()
            next_ref_meth = list()
            next_alt_meth = list()
        try:
            meth = 1 / (1. + math.exp(-float(row["log_lik_ratio"])))
        except OverflowError:
            meth = 0.0
        genotype = get_genotype(row["read_name"])
        if genotype == "ref":
            ref_meth.append(meth)
            ref_reads.add(row["read_name"])
        elif genotype == "alt":
            alt_meth.append(meth)
            alt_reads.add(row["read_name"])
        if row["start"] >= next_start:
            if genotype == "ref":
                next_ref_meth.append(meth)
                next_ref_reads.add(row["read_name"])
            elif genotype == "alt":
                next_alt_meth.append(meth)
                next_alt_reads.add(row["read_name"])

# print last result
if chr is not None:
    ref_mean = np.mean(ref_meth) if len(ref_meth) > 0 else "NA"
    ref_sd = np.std(ref_meth) if len(ref_meth) > 0 else "NA"
    alt_mean = np.mean(alt_meth) if len(alt_meth) > 0 else "NA"
    alt_sd = np.std(alt_meth) if len(alt_meth) > 0 else "NA"
    print("\t".join([chr, str(start), str(end), str(ref_mean), str(ref_sd), str(len(ref_meth)), str(len(
        ref_reads)), str(alt_mean), str(alt_sd), str(len(alt_meth)), str(len(alt_reads))]))
