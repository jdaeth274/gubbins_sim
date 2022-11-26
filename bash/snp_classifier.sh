#!/bin/bash 

set -e 

if [ $# -lt 4 ]
then
printf "This script needs four arguments not %s : \n" $#
echo "bash snp_classifier.sh <embl_dir> <gff_dir> <threads> <snp_classifier.R location> (OPTIONAL)<EMBL PREFIX LINE>"
echo "Both the directories should be just the embl_csv.csv files and the .gffs respectively"
echo ""

else

EMBL_DIR=$1
GFF_DIR=$2
THREADS=$3
SCRIPT_LOC=$4
#EMBL_PREFIX=$5
START=$SECONDS
ls -d "${EMBL_DIR}/"*.embl_csv.csv > embl_lists.txt
COUNTER=1
TOT_FILES=$(ls -d "${EMBL_DIR}/"*.embl_csv.csv | wc -l )
printf "Starting run through \n"

function my_func {
	BASIO=$(basename $1)
	MODEL_BASE=$(echo $BASIO | sed 's/\.embl_csv\.csv//g')
	GFF_FILE="${2}/${MODEL_BASE}.recombination_predictions.gff"
	if [ -f $GFF_FILE ]
	then
	Rscript --vanilla $3 --reccy-gff $GFF_FILE --branch-base $1 \
	--threads 1 --out "${MODEL_BASE}.classified_snps.csv" >> r_output_snp_class.txt 
	#printf "\r Finished on %s of %s embls" $4 $TOT_FILES
	COUNTER=$(( COUNTER + 1 ))
	else
		exit
	fi
}

export -f my_func
#parset line,GFF_DIR,SCRIPT_LOC,COUNTER 
START=$SECONDS
parallel -j $THREADS my_func ::: "${EMBL_DIR}/"*.embl_csv.csv ::: $GFF_DIR ::: $SCRIPT_LOC ::: $COUNTER

END=$(( SECONDS - START ))
printf "\n"
printf "Finished in %s (S) \n" $END
fi
