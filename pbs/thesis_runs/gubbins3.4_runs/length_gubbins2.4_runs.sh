#PBS -l select=1:ncpus=1:mem=64gb:cpu_type=rome
#PBS -l walltime=72:00:00
#PBS -J 1-100

set -e
module load anaconda3/personal

source activate gubbins2.4

## Switch to the run directory 
cd ~/../ephemeral/gubbins_testing/length_runs

## Get the run models 
MODEL=$( head -n $PBS_ARRAY_INDEX ~/gubbins_testing/length_tests.txt | tail -n 1)

SIMUL_LOC="/rds/general/user/jd2117/ephemeral/gubbins_testing/simul_data/length_tests/length_tests/"
REP=$( echo $MODEL | awk -F "-" '{print "rep_"$1}')
REP_DASH=$( echo $MODEL | awk -F "-" '{print "rep-"$1}')
SEQ_LENGTH=$(echo $MODEL | awk -F "-" '{print $2}')
FIRST_TREE="hybrid"
MAIN_TREE="hybrid"
MAR="hybrid"
MAIN_MODEL="${FIRST_TREE}-${MAIN_TREE}-${MAR}"
##Conditional statement for the reconstructions 
	cd "${REP}_runs"
	BASIO="sim-branch-0.1-rec-0.5-length-${SEQ_LENGTH}"
	ALNFILE="${SIMUL_LOC}${REP}_data/${BASIO}/${BASIO}-exact.aln"
	if [ -d $BASIO ]
	then
		cd $BASIO
	else
		mkdir -p $BASIO
		cd $BASIO
	fi

	## Get the alignment file

	if [ -d $MAIN_MODEL ]
	then
		rm -r $MAIN_MODEL
		mkdir $MAIN_MODEL
	else
		mkdir $MAIN_MODEL

	fi
	cd $MAIN_MODEL


	cp $ALNFILE ./
    TIME_PREFIX="${MAIN_MODEL}-${REP_DASH}-seq-${SEQ_LENGTH}"
	echo $PWD
		START_GUB=$SECONDS
		/usr/bin/time -v run_gubbins.py --prefix ${TIME_PREFIX} --use_time_stamp \
--threads 1 --verbose \
--tree_builder hybrid \
"${BASIO}-exact.aln" > gubbins_log 2>&1 
		END_GUB=$(( SECONDS - START_GUB ))
		PROC=$(cat /proc/cpuinfo | grep 'model name' | uniq)
		echo "Took this long for $TIME_PREFIX run, $END_GUB" > "../${MAIN_MODEL}_time.txt"
		echo "Used this processor for the runs: $PROC" >> "../${MAIN_MODEL}_time.txt"
		rm *.aln
echo "Finished the runs for this model -  ${MODEL}-${PBS_ARRAY_INDEX}" >> "./models_finished_old_gubb.txt"
