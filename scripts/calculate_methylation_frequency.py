#! /usr/bin/env python
# Modified version of script from nanopolish to aggregate metylation over genomic positions
# Adds -p option to aggregate by expected values (E(beta) =
# 1/(1+e^(-LLR))) instead of binary calls
from __future__ import print_function
import math
import sys
import csv
import argparse
from collections import namedtuple


def make_key(c, s, e):
    return c + ":" + str(s) + ":" + str(e)


class SiteStats:

    def __init__(self, g_size, g_seq):
        self.num_reads = 0
        self.posterior_methylated = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.group_size = g_size
        self.sequence = g_seq


def update_call_stats(key, num_called_cpg_sites, methylation, sequence):
    if key not in sites:
        sites[key] = SiteStats(num_called_cpg_sites, sequence)

    sites[key].num_reads += 1
    sites[key].called_sites += num_called_cpg_sites
    if methylation > 0:
        sites[key].called_sites_methylated += num_called_cpg_sites * methylation

parser = argparse.ArgumentParser(
    description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float,
                    required=False, default=2.5)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-s', '--split-groups', action='store_true')
parser.add_argument('-p', '--probabilistic', action='store_true')
parser.add_argument('--prior', type=float, default=0.5,
                    help="Prior probability of methylation")
args = parser.parse_args()
assert(args.call_threshold is not None)
if args.probabilistic:
    args.call_threshold = 0

sites = dict()

if args.input:
    in_fh = open(args.input)
else:
    in_fh = sys.stdin
csv_reader = csv.DictReader(in_fh, delimiter='\t')

for record in csv_reader:

    try:
        num_sites = int(record['num_motifs'])
    except KeyError:
        # backwards compatible
        num_sites = int(record['num_cpgs'])
    llr = float(record['log_lik_ratio'])

    # Skip ambiguous call
    if abs(llr) < args.call_threshold:
        continue
    sequence = record['sequence']

    if args.probabilistic:
        try:
            if llr > 700:
                llr = 700
            elif llr < -700:
                llr = -700
            methylation = 1 / (1 + (1 - args.prior) /
                               (args.prior * math.exp(llr)))
        except OverflowError:
            print(llr)
            raise
    else:
        methylation = 1 if llr > 0 else 0

    # if this is a multi-cpg group and split_groups is set, break up these
    # sites
    if args.split_groups and num_sites > 1:
        c = record['chromosome']
        s = int(record['start'])
        e = int(record['end'])

        # find the position of the first CG dinucleotide
        sequence = record['sequence']
        cg_pos = sequence.find("CG")
        first_cg_pos = cg_pos
        while cg_pos != -1:
            key = make_key(c, s + cg_pos - first_cg_pos,
                           s + cg_pos - first_cg_pos)
            update_call_stats(key, 1, methylation, "split-group")
            cg_pos = sequence.find("CG", cg_pos + 1)
    else:
        key = make_key(record['chromosome'], record['start'], record['end'])
        update_call_stats(key, num_sites, methylation, sequence)

# header
print("\t".join(["chromosome", "start", "end", "num_cpgs_in_group",
                 "called_sites", "called_sites_methylated",
                 "methylated_frequency", "group_sequence"]))

for key in sites:
    if sites[key].called_sites > 0:
        (c, s, e) = key.split(":")
        f = float(sites[key].called_sites_methylated) / sites[key].called_sites
        print("\t".join([str(x) for x in [c, s, e, sites[key].group_size,
                                          sites[key].called_sites,
                                          sites[key].called_sites_methylated,
                                          f, sites[key].sequence]]))
