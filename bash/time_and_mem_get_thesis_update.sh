#! /bin/bash/

#set -e 
## Point  this to a direcotry where the models are stored in the fasttree-iqtree-jar format, thats MODEL_DIR variable 
## Use the R script to get the memory usage and run time for each model. 

#GUBBINS_LOG=$1
OUT_FILE=$1
MODEL_FILE=$2
MODEL_DIR=$3
MODEL_NUM=1

if [ $# -ne 3 ]
then
echo "Need three arguments see usage below:"
echo ""
echo "bash time_and_mem_get.sh <out_name> <list_of_model_file> <results_dir>"
echo ""
exit
fi

if [ -f $OUT_FILE ]
then
echo "Output file already exists, overwriting this one completely"
rm $OUT_FILE
fi

cat $MODEL_FILE | while read line
do 
echo "$MODEL_NUM"
TIM_FILE="${MODEL_DIR}/${line}_time.txt"
echo $TIM_FILE
if [ -f $TIM_FILE ]
then 
cat $TIM_FILE >> $OUT_FILE
fi
CURRENT_GUBB="${MODEL_DIR}/${line}/gubbins_log"
echo $CURRENT_GUBB
grep "(wall clock)" $CURRENT_GUBB >> $OUT_FILE || echo "No time data for this run"
grep "Maximum resident set size" $CURRENT_GUBB >> $OUT_FILE || echo "No memory usage data for this run"
MODEL_NUM=$(( MODEL_NUM + 1 ))
done



