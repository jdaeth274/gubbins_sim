#PBS -l select=1:ncpus=1:mem=24gb
#PBS -l walltime=36:00:00
#PBS -J 1-180

set -e
module load anaconda3/personal

source activate gubbins_git

## Switch to the run directory 
cd ~/../ephemeral/gubbins_testing/length_runs

## Get the run models 
MODEL=$( head -n $PBS_ARRAY_INDEX ~/gubbins_testing/model_rep.txt | tail -n 1)

SIMUL_LOC="/rds/general/user/jd2117/ephemeral/gubbins_testing/simul_data/length_tests/"
REP=$( echo $MODEL | awk -F "-" '{print $4"_"$5"}')
REP_DASH=$( echo $MODEL | awk -F "-" '{print $4"-"$5"}'
FIRST_TREE=$( echo $MODEL | awk -F "-" '{print $1}')
MAIN_TREE=$( echo $MODEL | awk -F "-" '{print $2}')
MAR=$( echo $MODEL | awk -F "-" '{print $3}')
MAIN_MODEL="${FIRST_TREE}-${MAIN_TREE}-${MAR}"
##Conditional statement for the reconstructions 
cat ~/gubbins_testing/seq_lengths.txt | while read line 
do
	cd "${REP}_runs"
        BASIO="sim-branch-0.1-rec-0.5-length-${line}"
	ALNFILE="${SIMUL_LOC}${REP}_data/${BASIO}/${BASIO}.aln"
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
    TIME_PREFIX="${MAIN_MODEL}-${REP_DASH}-seq-${line}"

	if [ $MAR == "mar" ]
	then

		START_GUB=$SECONDS
		/usr/bin/time -v python ~/gubbins/python/run_gubbins.py --prefix ${TIME_PREFIX} --use-time-stamp --threads 1 --verbose \
		--tree-builder $MAIN_TREE --first-tree-builder $FIRST_TREE --mar \
		--model JC "${BASIO}.aln" > gubbins_log 2>&1 
		END_GUB=$(( SECONDS - START_GUB ))
		PROC=$(cat /proc/cpuinfo | grep 'model name' | uniq)
		echo "Took this long for $TIME_PREFIX run, $END_GUB" > "../${MAIN_MODEL}_time.txt"
		echo "Used this processor for the runs: $PROC" >> "../${MAIN_MODEL}_time.txt"
		rm *.aln
		cd ../../../


	else

		START_GUB=$SECONDS
		/usr/bin/time -v python ~/gubbins/python/run_gubbins.py --prefix ${TIME_PREFIX} --use-time-stamp --threads 1 --verbose \
		--tree-builder $MAIN_TREE --first-tree-builder $FIRST_TREE \
		--model JC "${BASIO}.aln" > gubbins_log 2>&1
		END_GUB=$(( SECONDS - START_GUB ))
		rm *.aln
		PROC=$(cat /proc/cpuinfo | grep 'model name' | uniq)
		echo "Took this long for ${TIME_PREFIX} run, $END_GUB" > "../${MAIN_MODEL}_time.txt"
		echo "Used this processor for the runs: $PROC" >> "../${MAIN_MODEL}_time.txt"
		cd ../../../
	fi
done

echo "Finished the runs for this model -  ${MODEL}-${PBS_ARRAY_INDEX}" >> "./models_finished_branch.txt"
