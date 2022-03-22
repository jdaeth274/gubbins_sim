#!/bin/bash

set -e
# Get the rep prefix
REP_NUM=$1


if [ -d "rep_${REP_NUM}_branchs" ]
then
	rm -r "rep_${REP_NUM}_branchs"
	mkdir "rep_${REP_NUM}_branchs"
else
	mkdir "rep_${REP_NUM}_branchs"
fi

if [ -d "rep_${REP_NUM}_trees" ]
then
	rm -r "rep_${REP_NUM}_trees"
	mkdir "rep_${REP_NUM}_trees"
else
	mkdir "rep_${REP_NUM}_trees"
fi

if [ -d "rep_${REP_NUM}_embls" ]
then	
	rm -r "rep_${REP_NUM}_embls"
	mkdir "rep_${REP_NUM}_embls"
else
	mkdir "rep_${REP_NUM}_embls"
fi

if [ -d "rep_${REP_NUM}_gffs" ]
then	
	rm -r "rep_${REP_NUM}_gffs"
	mkdir "rep_${REP_NUM}_gffs"
else
	mkdir "rep_${REP_NUM}_gffs"
fi


ls -d *-*-*/ > gubbins_list.txt

while read line <&3
do
cd $line
ls -d *-*-*-* > dir_conds.txt
while read sim <&4
do
cd $sim
cp *per_branch* ../../"rep_${REP_NUM}_branchs"
cp *node_labelled* ../../"rep_${REP_NUM}_trees"
cp *embl_csv* ../../"rep_${REP_NUM}_embls"
cp *recombination_predictions.gff* ../../"rep_${REP_NUM}_gffs"
cd ../
done 4< dir_conds.txt
echo $line
cd ../
done 3< gubbins_list.txt


