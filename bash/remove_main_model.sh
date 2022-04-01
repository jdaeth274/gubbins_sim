#!/bin/bash 

set -e 

if [ $# -ne 1 ]
then
printf "This needs a suffix, need the one extra argument \n"
printf "bash remove_main_model.sh <file_suffix> \n"
printf "\n"
else
ls sim-branch-*joint*$1 | while read line
do
	MODDY=${line:23}
	FIRST=${MODDY::1}
	if [ $FIRST == "-" ]
	then
		MODDY=${MODDY:1}
	fi
	echo $line
	echo $MODDY
	mv $line $MODDY
	printf "Done on: %s \n" $MODDY
done
fi

echo "Done"

