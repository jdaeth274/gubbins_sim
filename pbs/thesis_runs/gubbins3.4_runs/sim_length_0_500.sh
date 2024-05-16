#PBS -l select=1:ncpus=1:mem=6gb
#PBS -l walltime=24:00:00
#PBS -J 1-10

set -e
module load anaconda3/personal

source activate bioperl

cd ~/../ephemeral/gubbins_testing/simul_data/length_tests/length_tests

LENGTH=$(head -n $PBS_ARRAY_INDEX ~/gubbins_testing/seq_lengths.txt | tail -n 1)

#if [ $RECCERS == "0.1" ]
#then
#	echo "wassssssssup"
#else

for REP in {1..10}
do
cd "rep_${REP}_data"
## Get the simulation done
bash ~/gubbins_testing/gubbins_sim/bash/simul_script_poisson_JC_VARYING_LENGTH.sh \
0.1 0.5 \
~/gubbins_testing/simul_data/strep_53_uc.aln \
~/gubbins_testing/simul_data/26_phage.mfa \
$LENGTH

AIM_LENGTH=$(( LENGTH * 2 ))
NUM_ISOS=$(grep -c "^>" "sim-branch-0.1-rec-0.5-length-${LENGTH}.aln")
if [ $NUM_ISOS -ne $LENGTH ]
then 
printf "Rep %s Length %s had %s isolates \n" $REP $LENGTH $NUM_ISOS >> ../../altered_lengths.txt
fi
head -n $AIM_LENGTH "sim-branch-0.1-rec-0.5-length-${LENGTH}.aln" > "sim-branch-0.1-rec-0.5-length-${LENGTH}-exact.aln"



cd ../
## Make the IQtree
done

