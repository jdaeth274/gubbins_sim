#!/bin/bash

set -e 

## Running in the gubbins testing locales with all the rep_x_runs files 

ls -d rep_*_runs > rep_dirs.txt 

if [ -d jc_rep_10_csvs ]
then
	echo "Using already created csvs directory"
    sleep 3s
else
	mkdir jc_rep_10_csvs
fi

if [ -d jc_rep_10_trees ]
then
	echo "Using already created trees directory"
    sleep 3s
else
	mkdir jc_rep_10_trees
fi

if [ -d jc_rep_10_embls ]
then
	echo "Using already created embls directory"
    sleep 3s
else
	mkdir jc_rep_10_embls
fi

if [ -d jc_rep_10_gffs ]
then
	echo "Using already created gffs directory"
    sleep 3s
else
	mkdir jc_rep_10_gffs
fi


while read line <&3
do
    cd $line
    CURRENT_REP=$(echo $line | awk -F "_" '{print $2}')
    printf "Working on rep %s " $CURRENT_REP
    if [ -d "rep_${CURRENT_REP}_branchs" ]
    then
        cp "rep_${CURRENT_REP}_branchs"/* ../jc_rep_10_csvs
        printf "."
    else
        echo ""
        echo "No branch directory for rep $CURRENT_REP " 
    fi

    if [ -d "rep_${CURRENT_REP}_embls" ]
    then
        cp "rep_${CURRENT_REP}_embls"/* ../jc_rep_10_embls
        printf "."
    else
        echo "No embl directory for rep $CURRENT_REP "
    fi

    if [ -d "rep_${CURRENT_REP}_gffs" ]
    then
        cp "rep_${CURRENT_REP}_gffs"/* ../jc_rep_10_gffs
        printf "."
    else
        echo "No gff directory for rep $CURRENT_REP "
    fi

        if [ -d "rep_${CURRENT_REP}_trees" ]
    then
        cp "rep_${CURRENT_REP}_trees"/* ../jc_rep_10_trees
        printf ". Done \n"
    else
        echo "No tree directory for rep $CURRENT_REP "
    fi

    cd ../
    echo "Finished on $CURRENT_REP "
done 3< rep_dirs.txt



