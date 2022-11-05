#! /usr/local/bin/bash

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"
PARENT="$(dirname "$BASE_DIR")"
PERLDIR="${PARENT}/perl/"

BRANCH_RATE=$1
REC_RATE=$2
MIN_LENGTH=$3
DONOR_ALN=$4
PHAGE_ALN=$5

if [ $# -ne 5 ]
then
	echo "Incorrect number of argument need 5 have $# "
	echo "   "
	echo "bash <branch_rate> <rec_rate> <min_length> <donor_aln> <phage_aln>"
	echo ""
	exit 1
fi

dir="sim-branch-${BRANCH_RATE}-rec-${REC_RATE}"
if [ -d $dir ]
then
	rm -r $dir
	mkdir $dir
else
	mkdir $dir
fi
cd $dir
perl "${PERLDIR}generate_taxa_poisson_rec.pl" -a $DONOR_ALN \
-o $dir -n 100 -m $PHAGE_ALN -i 0 -b $BRANCH_RATE -r $REC_RATE \
-p $MIN_LENGTH
