for K in {1..10};
 do cd rep_10_data;
 ls -d /rds/general/user/jd2117/home/gubbins_testing/gubbins_sim/*/ > dir_list.txt;
 cat dir_list.txt | while read line;
 do
 cd $line 
 CONDY=$(basename $line)
 if [ -e orig-orig-orig*.tre ];
 then
 printf "${K}--${line}" >> ../../done_sims.txt;
 else
 printf "This dataset has no tree ${K}--${line}";
 fi;
 cd ../;
 done;
 cd ../;
 done 
