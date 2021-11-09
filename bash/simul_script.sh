#! /usr/local/bin/bash

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"
PARENT="$(dirname "$BASE_DIR")"
PERLDIR="${PARENT}/perl/"

BRANCH_RATE=$1
REC_RATE=$2
DONOR_ALN=$3
PHAGE_ALN=$4

if [ $# -ne 4 ]
then
	echo "Incorrect number of argument need 4 have $# "
	echo ""
	echo "bash <branch_rate> <rec_rate> <donor_aln> <phage_aln>
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
perl "${PERLDIR}generate_taxa.pl" -a $DONOR_ALN \
-o $dir -n 100 -m $PHAGE_ALN -i 0 -b $BRANCH_RATE -r $REC_RATE
#~nc3/Scripts/run_gubbins.py preliminary_run.aln

#VCF=`ls -1 *fixed_vers*vcf | tail -n 1 | tr -d "\n"`;

#~sh16/scripts/reportlabtest.py -M -t ${VCF%.vcf} -o preliminary_run.gubbins.pdf -q taxa ${VCF%.vcf}.tab
#~sh16/scripts/Find_recombinations_on_tree.py -a preliminary_run.aln -p preliminary_run.Find_rec -R -m 3 # This is the original gubbins script
#~sh16/scripts/reportlabtest.py -M -t preliminary_run.Find_rec_Final.tre -o preliminary_run.Find_rec.pdf -q taxa preliminary_run.Find_rec_rec.tab # A plotting of gubb results
#~nc3/Scripts/gubbins_evaluation.pl -t preliminary_run.summary  -s preliminary_run.Find_rec_SNPS_per_branch.tab -o test.out -g preliminary_run.Find_rec_rec.tab -a preliminary_run.aln # More summary statistics
