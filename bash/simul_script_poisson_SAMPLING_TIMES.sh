#! /usr/local/bin/bash

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"
PARENT="$(dirname "$BASE_DIR")"
PERLDIR="${PARENT}/perl/"

BRANCH_RATE=$1
REC_RATE=$2
DONOR_ALN=$3
PHAGE_ALN=$4
NUM_SEQS=$5
PARETO=$6
SAMPLE_RATE=$7

if [ $# -ne 7 ]
then
	echo "Incorrect number of argument need 5 have $# "
	echo "   "
	echo "bash <branch_rate> <rec_rate> <donor_aln> <phage_aln> <num_seqs> <minimum_recombination_length> <sample_rate>"
	echo ""
	exit 1
fi

dir="sim-branch-${BRANCH_RATE}-rec-${REC_RATE}-length-${NUM_SEQS}-pareto-${PARETO}"
if [ -d $dir ]
then
	rm -r $dir
	mkdir -p $dir
else
	mkdir -p $dir
fi
cd $dir
perl "${PERLDIR}generate_taxa_poisson_rec_SAMPLING.pl" -a $DONOR_ALN \
-o $dir -n $NUM_SEQS -m $PHAGE_ALN -i 0 -b $BRANCH_RATE -r $REC_RATE \
-p $PARETO -t $SAMPLE_RATE
