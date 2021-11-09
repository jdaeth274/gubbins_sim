#! /usr/local/bin/bash

dir="uppercase_runs";
mkdir ${dir};
cd ${dir};
perl ~/Dropbox/phd/gubbins_testing/analysis_code/generate_taxa.pl -a ../strep_alignment/split_kmers/strep_53_uc.aln \
-o 100_iso_run -n 100 -m ../26_phage.mfa -i 0 -b 0.1 -r 0.1
#~nc3/Scripts/run_gubbins.py preliminary_run.aln

#VCF=`ls -1 *fixed_vers*vcf | tail -n 1 | tr -d "\n"`;

#~sh16/scripts/reportlabtest.py -M -t ${VCF%.vcf} -o preliminary_run.gubbins.pdf -q taxa ${VCF%.vcf}.tab
#~sh16/scripts/Find_recombinations_on_tree.py -a preliminary_run.aln -p preliminary_run.Find_rec -R -m 3 # This is the original gubbins script
#~sh16/scripts/reportlabtest.py -M -t preliminary_run.Find_rec_Final.tre -o preliminary_run.Find_rec.pdf -q taxa preliminary_run.Find_rec_rec.tab # A plotting of gubb results
#~nc3/Scripts/gubbins_evaluation.pl -t preliminary_run.summary  -s preliminary_run.Find_rec_SNPS_per_branch.tab -o test.out -g preliminary_run.Find_rec_rec.tab -a preliminary_run.aln # More summary statistics
