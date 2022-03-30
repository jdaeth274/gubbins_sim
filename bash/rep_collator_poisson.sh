#!/bin/bash

set -e

# To be run within the gubbins_testing directory where the rep_x_runs data are 

SKIP_INDIV_COLLATE=$1

if [ $SKIP_INDIV_COLLATE == "No" ]
then
for K in {1..10}
do
	cd "rep_${K}_runs"
	bash ~/gubbins_testing/gubbins_sim/bash/rep_mover.sh $K /rds/general/user/jd2117/ephemeral/gubbins_testing/simul_data/poisson_data/"rep_${K}_data"
	cd ../
done
else
	echo "Skipping individual collation"
fi

bash ~/gubbins_testing/gubbins_sim/bash/rep_collator.sh
tar -czvf  jc_rep_10_csvs.tar.gz jc_rep_10_csvs/
cp jc_rep_10_csvs.tar.gz ~/
tar -czvf  jc_rep_10_gffs.tar.gz jc_rep_10_gffs/
cp jc_rep_10_gffs.tar.gz ~/
tar -czvf  jc_rep_10_embls.tar.gz jc_rep_10_embls/
cp jc_rep_10_embls.tar.gz ~/
tar -czvf  jc_rep_10_trees.tar.gz jc_rep_10_trees/
cp jc_rep_10_trees.tar.gz ~/

echo "Done !"

