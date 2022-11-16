#! /bin/bash

set -e 

cd ~/../ephemeral/gubbins_testing/pareto_20k_runs

mkdir -p pareto_20k_gubbins_trees
mkdir -p pareto_20k_gubbins_csvs



for K in {1..10}
do
    cd "rep_${K}_runs"
    cat dir_list.txt | while read line;
    do 
        cd $line
        while read model <&3
        do  
            cd $model 
            SIM=$(basename $line)
	    echo $PWD
            cp *per_branch* "../../../pareto_20k_gubbins_csvs/${model}-rep-${K}-${SIM}.per_branch_statistics.csv"
            cp *node_labelled* "../../../pareto_20k_gubbins_trees/${model}-rep-${K}-${SIM}.node_labelled.final_tree.tre" 
            printf "Done on ${K}-${SIM}-${model}\n"
            cd ../
        done 3< ~/gubbins_testing/gubbins_3.3_simple_complex.txt
        cd ../
    done 
    cd ../
done 

echo "Done on all isos!!!"
        
