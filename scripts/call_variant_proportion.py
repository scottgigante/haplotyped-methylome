#!/usr/bin/env python
# Calls haplotypes on each read based on results from nanopolish phase-reads as well as the original alignment
# Usage: ./call_variant_proportion.py -b /path/to/original.bam -p
# /path/to/phased.bam -v /path/to/variants.vcf

from __future__ import print_function
import sys
import pysam
import csv
import pickle
import functools
import multiprocessing
import argparse
import math

__min_coverage = 5
__override_ratio = 3

# fit exponential to log(1-correctness)
__slope = -0.1203
__intercept = -0.6927


def parse_args():
    parser = argparse.ArgumentParser(
        description="Phase nanopore reads from nanopolish output")
    parser.add_argument('-b', '--bam', required=True,
                        help="Path to the original bam file")
    parser.add_argument('-p', '--phased', required=True,
                        help="Path to the phased bam file from nanopolish phase-reads")
    parser.add_argument('-v', '--vcf', required=True,
                        nargs='+', help="Path to the variants file(s)")
    parser.add_argument('-d', '--detailed-output', action='store_true',
                        default=False, help="Print per-read, per-site calls to stdout")
    parser.add_argument('-o', '--outfile',
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-t', '--threads', type=int,
                        default=multiprocessing.cpu_count())
    args = parser.parse_args()
    return args


def get_read_id(read, contig, position, suppdb):
    if read.query_name in suppdb:
        for i, alignment in enumerate(suppdb[read.query_name]):
            chr, start, end = alignment
            if contig == chr and position >= start and position <= end:
                return "{}_{}".format(read.query_name, i)
        raise Exception("Position {}:{} not internal to any alignment of read {}".format(
            contig, position, read.query_name))
    else:
        return read.query_name


def check_column_support(vcf_rows, vcf_pos, column, reads, suppdb, type, detailed_output):
    if column is not None:
        valid_allele = [True if vcf_row['POS'] == vcf_pos and len(
            vcf_row['REF']) == 1 else False for vcf_row in vcf_rows]
        if column.reference_pos == vcf_pos and True in valid_allele and 'PASS' in [vcf_row['FILTER'] for vcf_row in vcf_rows]:
            for read in column.pileups:
                if not read.is_del:
                    read_name = read.alignment.query_name
                    read_pos = read.query_position
                    base = read.alignment.query_sequence[read_pos]
                    call = []
                    try:
                        qual = read.alignment.query_qualities[read_pos]
                        if type == "signal":
                            # nanopolish qualities scores should be offset by
                            # 35 - they never give anything below this
                            qual = max(0, qual - 35)
                        error = math.exp(__intercept + __slope * qual)
                    except TypeError:
                        if read.alignment.query_qualities is None:
                            # quality scores missing
                            error = 0
                        else:
                            # unknown error
                            raise
                    try:
                        reads[read_name]
                    except KeyError:
                        reads[read_name] = [0] * (len(vcf_rows) + 2)

                    # only have to check reference once
                    if base == vcf_rows[0]["REF"]:
                        reads[read_name][0] += 1 - error
                        call.append("ref")
                    else:
                        reads[read_name][0] += error
                    for i in range(len(vcf_rows)):
                        if vcf_rows[i]['FILTER'] == 'PASS':
                            match = False
                            if valid_allele[i]:
                                if base in vcf_rows[i]["ALT"]:
                                    match = True
                            elif "ref" in call and base == vcf_rows[0]["REF"]:
                                # no snp here, same as reference
                                match = True
                            if match:
                                reads[read_name][i + 1] += 1 - error
                                call.append("alt{}".format(i) if len(
                                    vcf_rows) > 1 else "alt")
                            else:
                                reads[read_name][i + 1] += error
                        else:
                            # treat as heterozygous
                            if valid_allele[i] and base in vcf_rows[i]['REF'] + vcf_rows[i]['ALT']:
                                reads[read_name][i + 1] += (1 - error) / 2
                            else:
                                reads[read_name][i + 1] += error
                    # add 1 to coverage
                    reads[read_name][-1] += 1
                    if len(call) > 0 and detailed_output:
                        qual = read.alignment.query_qualities[read_pos]
                        print(vcf_rows)
                        print(base)
                        print(error)
                        print(reads[read_name])
                        print("{}\t{}\t{}\t{}\t{}\t{}".format(
                            column.reference_name, column.reference_pos, read_name, type, ",".join(call), qual))
    return reads


def ref_support(reads):
    return set([read for read, counts in reads.items() if counts[0] > counts[1]])


def summarize_reads(reads):
    if len(reads) == 0:
        return None
    support = len(ref_support(reads))
    return float(support) / len(reads)


def load_vcf(handle):
    header = next(handle)
    while header.startswith("##"):
        header = next(handle)
    header = header.strip("#").strip("\n").split("\t")
    reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
    return reader


def pileup_step(ref_pos, original_pileup, phased_pileup):
    try:
        original_column = next(original_pileup)
        while original_column.reference_pos < ref_pos:
            original_column = next(original_pileup)
    except StopIteration:
        original_column = None
    try:
        phased_column = next(phased_pileup)
        while phased_column.reference_pos < ref_pos:
            phased_column = next(phased_pileup)
    except StopIteration:
        phased_column = None
    return original_column, phased_column


def process_row(vcf_rows, vcf_pos, original_pileup, original_reads, phased_pileup, phased_reads, suppdb, detailed_output):
    if 'PASS' in [vcf_row['FILTER'] for vcf_row in vcf_rows]:
        # process this snp
        original_column, phased_column = pileup_step(
            vcf_pos, original_pileup, phased_pileup)
        if not (original_column or phased_column):
            return None, original_reads, None, phased_reads
        else:
            original_reads = check_column_support(
                vcf_rows, vcf_pos, original_column, original_reads, suppdb, "base", detailed_output)
            phased_reads = check_column_support(
                vcf_rows, vcf_pos, phased_column, phased_reads, suppdb, "signal", detailed_output)
    return original_pileup, original_reads, phased_pileup, phased_reads


def call_contig(contig, vcf_fns, original_fn, phased_fn, suppdb, detailed_output):
    with pysam.AlignmentFile(original_fn, 'r') as original_bam, pysam.AlignmentFile(phased_fn, 'r') as phased_bam:
        vcf_handles = []
        vcf_readers = []
        for vcf_fn in vcf_fns:
            handle = open(vcf_fn, 'r')
            vcf_handles.append(handle)
            vcf_readers.append(load_vcf(handle))

        original_pileup = original_bam.pileup(contig)
        original_reads = dict()

        phased_pileup = phased_bam.pileup(contig)
        phased_reads = dict()

        if original_pileup or phased_pileup:
            # find right location in vcf - maybe an index is better for this
            vcf_rows = [next(vcf_reader) for vcf_reader in vcf_readers]
            for i, reader in enumerate(vcf_readers):
                if not vcf_rows[i]['CHROM'] == contig:
                    for vcf_row in reader:
                        if vcf_row['CHROM'] != contig:
                            # not there yet
                            continue
                        vcf_rows[i] = vcf_row
                        break
                vcf_rows[i]['POS'] = int(vcf_rows[i]['POS']) - 1
            # process the files
            while True:
                try:
                    vcf_pos = min(
                        [vcf_row['POS'] for vcf_row in vcf_rows if vcf_row['CHROM'] == contig])
                except ValueError:
                    # empty list, all on wrong contig
                    break
                original_pileup, original_reads, phased_pileup, phased_reads = process_row(
                    vcf_rows, vcf_pos, original_pileup, original_reads, phased_pileup, phased_reads, suppdb, detailed_output)
                if not (original_pileup or phased_pileup):
                    # out of reads
                    break
                num_failed = 0
                for i in range(len(vcf_rows)):
                    try:
                        if vcf_rows[i]['CHROM'] == contig and vcf_rows[i]['POS'] <= vcf_pos:
                            vcf_rows[i] = next(vcf_readers[i])
                            if vcf_rows[i]['CHROM'] != contig:
                                vcf_rows[i]['POS'] = float('inf')
                                num_failed += 1
                            else:
                                vcf_rows[i]['POS'] = int(
                                    vcf_rows[i]['POS']) - 1
                    except StopIteration:
                        vcf_rows[i]['CHROM'] = None
                        vcf_rows[i]['POS'] = float('inf')
                        num_failed += 1
                if num_failed == len(vcf_rows):
                    # no more variants
                    break
    return original_reads, phased_reads

args = parse_args()

with open("{}.suppdb".format(args.bam), 'rb') as handle:
    suppdb = pickle.load(handle)

func = functools.partial(call_contig, vcf_fns=args.vcf, original_fn=args.bam,
                         phased_fn=args.phased, suppdb=suppdb, detailed_output=args.detailed_output)
with pysam.AlignmentFile(args.bam, 'r') as bam:
    original_reads, phased_reads = dict(), dict()
    contigs = [bam.get_reference_name(i) for i in range(bam.nreferences)]
if args.detailed_output:
    print("\t".join(["chr", "pos", "read_name", "caller", "genotype", "qual"]))

threads = 1 if args.detailed_output else min(args.threads, len(contigs))
if threads > 1:
    pool = multiprocessing.Pool(threads)
    data = pool.imap_unordered(func, contigs)
    pool.close()
    pool.join()
else:
    data = [func(contig) for contig in contigs]

for original_contig_reads, phased_contig_reads in data:
    original_reads.update(original_contig_reads)
    phased_reads.update(phased_contig_reads)

reads = set(original_reads.keys()).union(set(phased_reads.keys()))
signal_alt_header = "\t".join(["signal_alt{}".format(
    i + 1) for i in range(len(args.vcf))]) if len(args.vcf) > 1 else "signal_alt"
base_alt_header = "\t".join(["base_alt{}".format(
    i + 1) for i in range(len(args.vcf))]) if len(args.vcf) > 1 else "base_alt"
print("\t".join(['read', 'genotype', 'signal_ref', signal_alt_header,
                 'signal_coverage', 'base_ref', base_alt_header, 'base_coverage']), file=args.outfile)
for read in reads:
    try:
        signal_calls = phased_reads[read]
    except KeyError:
        signal_calls = [0] * (len(args.vcf) + 2)
    try:
        base_calls = original_reads[read]
    except KeyError:
        base_calls = [0] * (len(args.vcf) + 2)

    base_coverage = base_calls[-1]
    base_calls = base_calls[:-1]
    signal_coverage = signal_calls[-1]
    signal_calls = signal_calls[:-1]
    base_max = max(base_calls)
    signal_max = max(signal_calls)

    if signal_coverage == 0 or signal_calls.count(signal_max) > 1:
        signal_call = "fail"
        signal_ratio = 0
    else:
        signal_call = signal_calls.index(signal_max)
        signal_ratio = signal_max / signal_coverage
        if signal_call == 0:
            signal_call = "ref"
        else:
            signal_call = "alt{}".format(
                signal_call) if len(args.vcf) > 1 else "alt"

    if base_coverage == 0 or base_calls.count(max(base_calls)) > 1:
        base_call = "fail"
        base_ratio = 0
    else:
        base_call = base_calls.index(max(base_calls))
        base_ratio = base_max / base_coverage
        if base_call == 0:
            base_call = "ref"
        else:
            base_call = "alt{}".format(base_call) if len(
                args.vcf) > 1 else "alt"

    if base_coverage < __min_coverage and signal_coverage < __min_coverage:
        genotype = "fail"
    elif base_call == signal_call:
        genotype = base_call
    elif signal_call == "fail":
        genotype = base_call
    elif base_call == "fail":
        genotype = signal_call
    elif base_coverage > __override_ratio * signal_coverage and base_call != "fail":
        genotype = base_call
    elif signal_coverage > __override_ratio * base_coverage and signal_call != "fail":
        genotype = signal_call
    elif base_ratio - 0.5 > __override_ratio * (signal_ratio - 0.5) and base_ratio > 0.5:
        genotype = base_call
    elif signal_ratio - 0.5 > __override_ratio * (base_ratio - 0.5) and signal_ratio > 0.5:
        genotype = signal_call
    else:
        genotype = "fail"

    if read in suppdb:
        for j in range(len(suppdb[read])):
            supp_read = "{}_{}".format(read, j)
            print("\t".join([supp_read, genotype, "\t".join([str(i) for i in signal_calls]), str(
                signal_coverage), "\t".join([str(i) for i in base_calls]), str(base_coverage)]), file=args.outfile)
    else:
        print("\t".join([read, genotype, "\t".join([str(i) for i in signal_calls]), str(
            signal_coverage), "\t".join([str(i) for i in base_calls]), str(base_coverage)]), file=args.outfile)
