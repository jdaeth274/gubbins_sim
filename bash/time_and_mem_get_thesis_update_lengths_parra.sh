#! /bin/bash/

#set -e 
## Point  this to the parent directory of the rep runs  
## Use the R script to get the memory usage and run time for each model. 

#GUBBINS_LOG=$1
OUT_FILE=$1
MODEL_FILE=$2
MODEL_NUM=1

if [ $# -ne 2 ]
then
echo "Need two arguments see usage below:"
echo ""
echo "bash time_and_mem_get.sh <out_name> <list_of_model_file> "
echo ""
exit
fi

if [ -f $OUT_FILE ]
then
echo "Output file already exists, overwriting this one completely"
rm $OUT_FILE
fi

##Loop through the rep runs 
for K in {1..10};
do 
cd "rep_${K}_runs"
## Now we need to go through each of the sim conditions 
cd "sim-branch-0.1-rec-0.5-length-500"
cat $MODEL_FILE | while read line
do 
echo "$MODEL_NUM"
TIM_FILE="${line}_time.txt"
echo $TIM_FILE
if [ -f $TIM_FILE ]
then 
cat $TIM_FILE >> $OUT_FILE
CURRENT_GUBB="${line}/gubbins_log"
echo $CURRENT_GUBB
grep "(wall clock)" $CURRENT_GUBB >> $OUT_FILE || echo "No time data for this run"
grep "Maximum resident set size" $CURRENT_GUBB >> $OUT_FILE || echo "No memory usage data for this run"
fi
MODEL_NUM=$(( MODEL_NUM + 1 ))
done
cd ../../
done 


