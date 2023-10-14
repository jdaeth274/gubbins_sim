#PBS -l select=1:ncpus=1:mem=32gb:cpu_type=rome
#PBS -l walltime=72:00:00
#PBS -J 1-240

set -e
module load anaconda3/personal

source activate gubbins_git

## Switch to the run directory 
cd ~/../ephemeral/gubbins_testing/time_runs

## Get the run models 
MODEL=$( head -n $PBS_ARRAY_INDEX ~/gubbins_testing/gubbins_3.3_models.txt | tail -n 1)

SIMUL_LOC="/rds/general/user/jd2117/ephemeral/gubbins_testing/simul_data/time_tests/"
REP=$( echo $MODEL | awk -F "-" '{print $4"_"$5}')
REP_DASH=$( echo $MODEL | awk -F "-" '{print $4"-"$5}')
FIRST_TREE=$( echo $MODEL | awk -F "-" '{print $1}')
MAIN_TREE=$( echo $MODEL | awk -F "-" '{print $2}')
MAR=$( echo $MODEL | awk -F "-" '{print $3}')
MAIN_MODEL="${FIRST_TREE}-${MAIN_TREE}-${MAR}"
##Conditional statement for the reconstructions 
for K in {1..2}; 
do

if [ $K -eq 1 ]
then 
    PARETO="0"
    REC="1"
else
    PARETO="20000"
    REC="0.1"
fi

cd "${REP}_runs"
BASIO="sim-branch-0.1-rec-${REC}-length-100-pareto-${PARETO}"
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
TIME_PREFIX="${MAIN_MODEL}-${REP_DASH}-seq-100-pareto-${PARETO}"
echo $PWD
if [ $MAR == "mar" ]
then

	START_GUB=$SECONDS
    FIRST_TREE=$(echo $FIRST_TREE | sed 's/_/-/g')
    /usr/bin/time -v python ~/gubbins/python/run_gubbins.py --prefix ${TIME_PREFIX} --use-time-stamp --threads 1 --verbose \
	--tree-builder $MAIN_TREE --first-tree-builder $FIRST_TREE --mar \
	--model JC "${BASIO}-exact.aln" > gubbins_log 2>&1 
	END_GUB=$(( SECONDS - START_GUB ))
	rm *.aln
	cd ../../../

else

	START_GUB=$SECONDS
    FIRST_TREE=$(echo $FIRST_TREE | sed 's/_/-/g')

	/usr/bin/time -v python ~/gubbins/python/run_gubbins.py --prefix ${TIME_PREFIX} --use-time-stamp --threads 1 --verbose \
	--tree-builder $MAIN_TREE --first-tree-builder $FIRST_TREE \
	--model JC "${BASIO}-exact.aln" > gubbins_log 2>&1
	END_GUB=$(( SECONDS - START_GUB ))
	rm *.aln
	cd ../../../
fi

done 
echo "Finished the runs for this model -  ${MODEL}-${PBS_ARRAY_INDEX}" >> "./models_finished_branch.txt"
