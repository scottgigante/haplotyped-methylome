# Usage: python mask_genome_variants.py reference.fa reference.masked.fa
# variants.vcf
from __future__ import print_function
from Bio import SeqIO
from Bio import Seq
import re
import sys
import csv
import argparse


def load_vcf(handle):
    header = next(handle)
    while header.startswith("##"):
        header = next(handle)
    header = header.strip("#").strip("\n").split("\t")
    reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
    return reader


__ambiguity_codes = {
    'N': {'A', 'C', 'G', 'T'},
    'V': {'A', 'C', 'G'},
    'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'},
    'B': {'C', 'G', 'T'},
    'M': {'A', 'C'},
    'K': {'G', 'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'W': {'A', 'T'},
    'S': {'C', 'G'},
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'}
}


def get_ambiguity_base(bases):
    bases = set(bases)
    for ambiguity_base, possible_bases in __ambiguity_codes.items():
        if bases == possible_bases:
            return ambiguity_base
    raise Exception("Unrecognises bases: ['" + "','".join(list(bases)) + "']")


def disambiguate_base(base):
    base = str(base)
    return list(__ambiguity_codes[base])


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
    region = ''.join(region.split())  # remove whitespace
    region = re.split(':|-', region)
    start = 0
    end = None
    if len(region) < 1:
        raise argparse.ArgumentTypeError(
            "Region must specify a reference name")
    elif len(region) > 3:
        raise argparse.ArgumentTypeError("Region format not recognized")
    else:
        contig = region[0]
        try:
            start = int(re.sub(",|\.", "", region[1]))
        except IndexError:
            pass
        try:
            end = int(re.sub(",|\.", "", region[2]))
        except IndexError:
            pass
        finally:
            return contig, start, end


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='filename of input genome fasta')
    parser.add_argument('-o', '--output', required=True,
                        help='filename of output masked genome fasta')
    parser.add_argument('-v', '--vcf', required=True,
                        help='filename of vcf listing snps')
    parser.add_argument('-a', '--alternate-only', default=False, action='store_true',
                        help='replace reference base with alternate allele, rather than retaining both (Default: false)')
    parser.add_argument('-u', '--update', default=False, action='store_true',
                        help='Update an already masked genome with another VCF file (Default: false)')
    parser.add_argument('-r', '--region', type=parse_region, nargs="+", metavar="CHR:START-END",
                        help="Name and range of contigs to mask (Default: entire genome)")
    parser.add_argument('-b', '--boundary', type=int, default=0, metavar="INT",
                        help="Width of region boundary where both alleles are included")
    return parser.parse_args()


def contig_to_int(contig):
    try:
        return(int(contig))
    except ValueError:
        if contig == 'X':
            return max_autosome + 1
        elif contig == 'Y':
            return max_autosome + 2
        elif contig == 'MT':
            return max_autosome + 3
        else:
            return max_autosome + 4


def cmp_contigs(contig):
    return(contig_to_int(contig), contig)


args = parse_args()
regions = dict()
# format: { chr : [start, end, alternate_only] }
if args.region:
    for region in args.region:
        try:
            regions[region[0]]
        except KeyError:
            regions[region[0]] = []
        region_start, region_end = region[1], region[2]
        if region_start > 0:
            region_start = region_start + args.boundary // 2
        if region_end is not None:
            region_end = region_end - args.boundary // 2
        regions[region[0]].append(
            [region_start, region_end, args.alternate_only])
        if args.boundary > 1:
            if region[1] > 0:
                regions[region[0]].append(
                    [region[1] - args.boundary // 2, region[1] + args.boundary // 2, False])
            if region[2] is not None:
                regions[region[0]].append(
                    [region[2] - args.boundary // 2, region[2] + args.boundary // 2, False])

genome = dict()
with open(args.input, 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record.seq = record.seq.tomutable()
        genome[record.id] = record

with open(args.vcf, 'r') as handle:
    vcf = load_vcf(handle)
    for row in vcf:
        snp = True
        ref_base = [row['REF']]
        alt_bases = row['ALT'].split(',')
        for base in alt_bases + ref_base:
            if len(base) > 1:
                snp = False
        if not snp:
            # don't touch indels
            continue

        pos = int(row['POS']) - 1
        if regions:
            alternate_only = None
            try:
                for r in regions[row['CHROM']]:
                    if pos >= r[0] and pos < r[1]:
                        alternate_only = r[2]
            except KeyError:
                pass
            if alternate_only is None:
                # no region hit
                continue
        else:
            # genome wide
            alternate_only = args.alternate_only

        record = genome[row['CHROM']]

        bases = alt_bases
        if not alternate_only:
            if args.update:
                bases = bases + disambiguate_base(record.seq[pos])
            else:
                bases = bases + ref_base

        record.seq[pos] = get_ambiguity_base(bases)
        genome[row['CHROM']] = record

contigs = list()
max_autosome = 0
for contig in genome.keys():
    contigs.append(contig)
    try:
        if int(contig) > max_autosome:
            max_autosome = int(contig)
    except ValueError:
        pass

contigs = sorted(contigs, key=cmp_contigs)

with open(args.output, 'w') as handle:
    for contig in contigs:
        SeqIO.write(genome[contig], handle, "fasta")
