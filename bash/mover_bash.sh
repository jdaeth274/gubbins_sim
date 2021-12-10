#!/bin/bash

set -e


if [ -d rec_csvs ]
then
	rm -r rec_csvs
	mkdir rec_csvs
else
	mkdir rec_csvs
fi

if [ -d node_trees ]
then
	rm -r node_trees
	mkdir node_trees
else
	mkdir node_trees
fi

ls -d *-*-*/ > gubbins_list.txt

while read line <&3
do
cd $line
ls -d *branch* > dir_conds.txt
while read sim <&4
do
cd $sim
cp *per_branch* ../../rec_csvs
cp *node_labelled* ../../node_trees
cd ../
done 4< dir_conds.txt
echo $line
cd ../
done 3< gubbins_list.txt


