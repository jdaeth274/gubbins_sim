#! /bin/bash

set -e 

cd ~/../ephemeral/gubbins_testing/pareto_0_runs

mkdir -p pareto_0_gubbins_trees
mkdir -p pareto_0_gubbins_csvs



for K in {1..10}
do
    cd "rep_${K}_runs"
    cat dir_list.txt | while read line;
    do 
        cd $line
        while read model <&3
        do  
            cd $model 
            cp *per_branch* ../../../pareto_0_gubbins_csvs
            cp *node_labelled* ../../pareto_0_gubbins_trees
            printf "Done on ${K}-${line}-${model}\n"
        done 3< ~/gubbins_testing/gubbins_3.3_simple_complex.txt
    done 
done 

echo "Done on all isos!!!"
        
