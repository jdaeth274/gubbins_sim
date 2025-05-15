#! /bin/bash

set -e 

cd ~/../ephemeral/gubbins_testing/pareto_0_runs

mkdir -p pareto_0_gubbins_trees
mkdir -p pareto_0_gubbins_csvs
mkdir -p pareto_0_gubbins_gff



for K in {1..10}
do
    cd "rep_${K}_runs"
    ls -d * > dir_list.txt
    cat dir_list.txt | while read line;
    do 
        cd $line
        while read model <&3
        do  
            cd $model 
            SIM=$(basename $line)
            if [ -e "${model}-rep-${K}.final_tree.tre"]
            then
                cp *per_branch* "../../../pareto_0_gubbins_csvs/${model}-rep-${K}-${SIM}-JC-JC.per_branch_statistics.csv"
                cp *node_labelled* "../../../pareto_0_gubbins_trees/${model}-rep-${K}-${SIM}-JC-JC.node_labelled.final_tree.tre" 
                cp *recombination_predictions.gff "../../../pareto_0_gubbins_gff/${model}-rep-${K}-${SIM}-JC-JC.recombination_predictions.gff"
                printf "Done on ${K}-${SIM}-${model}\n"
            else
                printf "Error on ${K}-${SIM}-${model}\n"
                printf "${K}-${SIM}-${model}\n" >> ../../../error_runs.txt
            cd ../
        done 3< ~/gubbins_testing/gubbins_3.3_simple_complex.txt
        cd ../
    done 
    cd ../
done 

echo "Done on all isos!!!"
        
