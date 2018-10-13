set -x

INFILE=$1
FASTQ=$2
BAM=$3
VCF=$4
TMP=$INFILE.tmp
mkdir -p $TMP

PHASED="${FASTQ}.phased.sorted.bam"
METH="${FASTQ}.methylation.tsv"

samtools index $BAM &
samtools index $PHASED &
wait

# build suppdb
python build_supplementary_db.py $BAM --summary ${INFILE}.summary.tsv &
python build_supplementary_db.py $PHASED &
wait

# split methylation and calculate summaries
METH_SPLIT=${FASTQ}.methylation.split_supplementary.tsv
PHASED_SUMMARY=${FASTQ}.phased.tsv
python split_methylation_by_alignment.py $BAM $METH > $METH_SPLIT & echo $! > ${TMP}/split.pid
python calculate_methylation_frequency.py -i $METH -p > ${FASTQ}.methylation.summary.tsv &
python call_variant_proportion.py -b $BAM -p $PHASED -v $VCF -o $PHASED_SUMMARY & echo $! > ${TMP}/phase.pid

# sort methylation
wait $(cat $TMP/split.pid)
METH_SITE=${FASTQ}.methylation.sorted.by_site.tsv
mkdir -p ~/tmp
tail -n +2 $METH_SPLIT | sort -k4,4 -k1,1 -k2n,2n | cat <(head -n 1 $METH) - > ${FASTQ}.methylation.sorted.by_read.tsv & echo $! > ${TMP}/read.pid
tail -n +2 $METH_SPLIT | sort -k1,1 -k2n,2n -k4,4 | cat <(head -n 1 $METH) - > $METH_SITE & echo $! > ${TMP}/site.pid

RDATA=../RData/$(basename $INFILE)/
wait $(cat ${TMP}/read.pid)
wait $(cat ${TMP}/phase.pid)
Rscript fit_reads.R $INFILE $RDATA &
wait $(cat ${TMP}/site.pid)
python split_methylation_by_haplotype.py -m $METH_SITE -p $PHASED_SUMMARY & echo $! > ${TMP}/split.pid

wait $(cat ${TMP}/split.pid)
\ & echo $! > ${TMP}/ref.pid
python calculate_methylation_frequency.py -i ${METH_SITE}.alt.tsv -p > ${FASTQ}.methylation.alt.summary.tsv & echo $! > ${TMP}/alt.pid

wait

rm -rf $TMP
set +x
