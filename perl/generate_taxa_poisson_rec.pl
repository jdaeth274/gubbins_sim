#! /software/perl-5.8.8/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long;
use Math::Random qw(random_poisson);


# parse and check command line input

my $help = 0;
my $aln_file;
my $mge_file;
my $outfile = "evolved_taxa";
my $n_taxa = 10;
my $branch_prob = 0.05;
my $rec_prob = 0.1;
my $ins_prob = 0.1;
my $debug = 0;

GetOptions (
	"h" => \$help,
	"a=s" => \$aln_file,
	"m=s" => \$mge_file,
	"o=s" => \$outfile,
	"n=s" => \$n_taxa,
	"b=s" => \$branch_prob,
	"r=s" => \$rec_prob,
	"i=s" => \$ins_prob,
	"v" => \$debug
);

if ($help == 1) {
	croak(); exit(0);
}

unless (defined($aln_file) && -e $aln_file && defined($mge_file) && -e $mge_file) {
	print STDERR "Need alignment file and MGE file!\n"; croak(); exit(1);
}

unless ($branch_prob >= 0 && $rec_prob >= 0 && $ins_prob >= 0) {
	print STDERR "branch, recombination and insertion probabilities must be positive rates in the range 0-100!\n"; croak(); exit(1);
}

# parse genome alignment

my $ancestor;
my $ancestor_name;
my %donor;

my $alnObject = Bio::AlignIO->new(-file => $aln_file, -format => "Fasta"); # Takes in a list of alignments
my $alignment = $alnObject->next_aln;

print STDERR "Parsing whole genome alignment...";

foreach my $seq ($alignment->each_seq) {  ## Reading in the alignment, first seq used as ancestor, others donors
	my $name = $seq->id;
	my $dna = $seq->seq;
	if (defined($ancestor) && length($ancestor) > 0) {
		$donor{$name} = $dna;
	} else {
		$ancestor_name = $name;
		$ancestor = $dna;
	}
}

print STDERR "\nGenerating alignment/sequence translation...";

# create translation which has a unique index for each position in the alignment

my %translate;
my %back_translate;
my $i = 1;
my $position = 1;
# Dictionary to translate and back translate due to MGE insertion
for (my $i = 1;$i <= length($ancestor);$i++) {$translate{$i} = $position; $position++;}

back_translate(); # Gets the original positions in the alignment

my $i_max = $position-1;

print STDERR "\nParsing MGE file...";

# parse MGE file

my %mge;

my $mgeObject = Bio::SeqIO->new(-file => $mge_file, -format => "Fasta");
# Separate fasta file of phages, when insert everything else gets gaps
while (my $seq = $mgeObject->next_seq) {
	my $name = $seq->id;
	my $dna = $seq->seq;
	$mge{$name} = $dna;
}

print STDERR "\nSimulating divergence…\n";

# run simulation of divergence

open LOG, "> $outfile.log" or print STDERR "Cannot open output file $outfile.log!\n" and croak() and exit(1);
open ALN, "> $outfile.aln" or print STDERR "Cannot open output file $outfile.aln!\n" and croak() and exit(1);
open SEQ, "> $outfile.mfa" or print STDERR "Cannot open output file $outfile.seq!\n" and croak() and exit(1);
open TAB, "> $outfile.tab" or print STDERR "Cannot open output file $outfile.tab!\n" and croak() and exit(1);
open ATAB, "> $outfile.aln.tab" or print STDERR "Cannot open output file $outfile.aln.tab!\n" and croak() and exit(1);
open TREE, "> $outfile.tree" or print STDERR "Cannot open output file $outfile.tree!\n" and croak() and exit(1);
open SUM, "> $outfile.summary" or print STDERR "Cannot open output file $outfile.summary!\n" and croak() and exit(1);

my %summary;
my %taxa;
my $index = 0;
$taxa{$index} = $ancestor;

my %event_log;
my $e_num = 0;
my %taxon_events;
my %event_pos;
my %event_type;
my %event_taxa;
my %log_insertions;
my %superinfection;	# do not allow a taxon to have two phage inserting at the same base
my %recsnp; # record whether it is possible to detect recombiantions or not

my %diverge;
my $treestring = "(taxon_0);";
my %exclude; # prevent recombinations overlapping point mutations or other recombinations on the same branch esp. in the terminal round of divergence
## Initiating the tree with the ancestor, first seq in aln
while (scalar(keys %taxa) < $n_taxa) {
	foreach my $taxon (keys %taxa) {
		$taxa{$taxon} = mutate($taxon);		# introduce a point mutation into every sequence at every timestep
		$diverge{$taxon}++;
		my $score = random_poisson(1, $rec_prob);
		if ($score > 0) {
			my $rec_counter = 0;
			while ($rec_counter < $score){
				$taxa{$taxon} = hom_rec($taxon);
				$rec_counter += 1;
			}										# introduce a homologous recombination. Function selects random start and length and replacing ancestor with donor
		}
		
		if ($ins_prob < 1) { # Phage/MGE insertions
			my $score = rand(1);
			if ($score <= $ins_prob) {
				phage_ins($taxon);		# introduce a phage into a sequence
			}
		} else {
			for (my $p = 1; $p <= $ins_prob**2; $p++) { # Again if > 1 square to get insertion prob
				my $score = rand($ins_prob**2);
				if ($score <= $ins_prob) {
					phage_ins($taxon);
				}
			
			}
		}
			
		my $score = rand(1);
		if ($score <= $branch_prob) { # Whether a new lineage is formed, sampling density, lower score more SNPs accumulate before detect branching
			$index++;
			$taxa{$index} = $taxa{$taxon};
			copy_event_log($taxon, $index); # Copies all the evo history into the new taxon
			my $to_find; my $to_replace;
			if ($index == 1) {
				$to_find = "taxon_"."$taxon";
				$to_replace  = "taxon_"."$taxon".",taxon_"."$index";
			} else {
				$to_find = "(taxon_"."$taxon)(,|"."\\)\;|"."\\"."):)"; # Updates the tree with a newick string
				$treestring =~ /$to_find/;
				if ($2 eq ',') {
					$to_replace  = "\(taxon_"."$taxon".",taxon_"."$index\):$diverge{$taxon}"."$2";
				} elsif ($2 eq ');') {
					$to_replace  = "\(taxon_"."$taxon".",taxon_"."$index\):$diverge{$taxon}"."$2";
				} elsif ($2 eq '):') {
					$to_replace  = "\(taxon_"."$taxon".",taxon_"."$index\):$diverge{$taxon}"."$2";
				}
			}
			undef(%{$exclude{$taxon}});	# mutations can now be distinguished
			$diverge{$taxon} = 0; # How many point mutations accumulating in taxon since ancestor
			$diverge{$index} = 0;
			$treestring =~ s/$to_find/$to_replace/g;
			print LOG "Taxa $index and $taxon diverge here\n";
		}
	}
}

# run one final round of divergence after the taxa threshold is reached; add in terminal branch lengths on tree

foreach my $taxon (keys %taxa) {
	$taxa{$taxon} = mutate($taxon);		# introduce a point mutation into every sequence at every timestep
	$diverge{$taxon}++;
	my $score = random_poisson(1, $rec_prob);
		if ($score > 0) {
			my $rec_counter = 0;
			while ($rec_counter < $score){
				$taxa{$taxon} = hom_rec($taxon);
				$rec_counter += 1;
			}										# introduce a homologous recombination. Function selects random start and length and replacing ancestor with donor
		}
		
	$score = rand(1);
	if ($score <= $ins_prob) {
		phage_ins($taxon);		# introduce a phage into a sequence
	}
	my $to_find = "(taxon_"."$taxon)(,|"."\\"."))";
	if ($treestring =~ /$to_find/) {
		my $to_replace = "$1"."\:"."$diverge{$taxon}"."$2";
		$treestring =~ s/$to_find/$to_replace/g;
	}
}

$treestring =~ s/,,/,/g;
$treestring =~ s/\,\)/\)/g;

print TREE "$treestring\n";

close TREE;

# print output multifasta;

print STDERR "Printing output!\n";

my $anc_bases = $ancestor; # INitial seq in the aln
$anc_bases =~ s/-|N|n//g;
print SEQ ">$ancestor_name\n$anc_bases\n";

my %gap_count;

foreach my $taxon (keys %taxa) {
	my $s; my $g = 0;
	my @s = split(//,$taxa{$taxon});
	foreach my $b (@s) {
		if ($b =~ /n|N|-/) {
			$gap_count{$g}++;
		} else {
			$s.=$b;
		}
		$g++;
	}
	print SEQ ">taxon_"."$taxon\n$s\n";
}

close SEQ; # Printing out the alignment

# establish relationship between taxon_0/ancestor and alignment

my %Otrans;
my $Opos = 0;
my $Oindex = 0;

my @t0 = split(//,$taxa{0});

foreach my $b (@t0) {
	if (!(defined($gap_count{$Oindex})) || $gap_count{$Oindex} < scalar(keys %taxa)) {
		$Opos++;
	}
	$Oindex++;
	$Otrans{$Oindex} = $Opos;
}
# Looking at positions in ancestor relative to the alignment, gap count incremented in evolution, where each pos in aln is relative to original ancestor
print STDERR "Relative to taxa alignment, maximum position $Opos maximum index $Oindex\n";

undef(@t0);

# print output alignment sans ancestor

foreach my $taxon (keys %taxa) {
	my $s; my $g = 0;
	my @s = split(//,$taxa{$taxon});	
	foreach my $b (@s) {
		if (!(defined($gap_count{$g})) || $gap_count{$g} < scalar(keys %taxa)) { # counting gap count in the column, if all gaps remove, could be MGE insert/delete
			$s.=$b;
		}
		$g++;
	}
	print ALN ">taxon_"."$taxon\n$s\n";
}

close ALN;

# establish relationship between reference sequence and alignment

my @anc = split(//,$ancestor);

my %rtrans;
my $rpos = 0;
my $rindex = 0;

while (@anc) {
	my $base = shift(@anc);
	unless ($base eq '-' || $base eq 'N' || $base eq 'n') {
		$rpos++;
	}
	$rindex++;
	$rtrans{$rindex} = $rpos;
}
# Once removed all gaps, shift the indices to make sense with the ancestor
# print debugging

if ($debug == 1) {
	open DEBUG, "> $outfile.debug";
	
	foreach my $p (sort {$translate{$a}<=>$translate{$b}} keys %translate) {
		print DEBUG "$p\t";
	}
	
	print DEBUG "\n";
	
	foreach my $p (sort {$translate{$a}<=>$translate{$b}} keys %translate) {
		print DEBUG "$translate{$p}\t";
	}
	
	print DEBUG "\n";
	
	foreach my $p (sort {$translate{$a}<=>$translate{$b}} keys %translate) {
		if (defined($rtrans{$p})) {
			print DEBUG "$rtrans{$p}\t";
		} else {
			print DEBUG "\t";
		}
	}
	
	print DEBUG "\n";
	
	foreach my $taxon (keys %taxa) {
		my @t = split(//,$taxa{$taxon});
		$" = "\t";
		
		print DEBUG ">$taxon\n@t\n";
	}
}

# end debug

# print summary and tab files

$" = ",";

print SUM "Event#\tType\tStart\tEnd\tTaxa\n";
# printing a log summary
foreach $e_num (keys %event_log) {
	my @extant;
	foreach my $t (@{$event_taxa{$e_num}}) {				# check if the mutation remains observable in any sequences
		if ($taxon_events{$t}{$e_num} == 1) {
			push(@extant,$t);
		}
	}
	if ($#extant > -1) {
		print SUM "$e_num\t$event_type{$e_num}"; # Print the event number then the event type
		if ($event_type{$e_num} eq 'S' || $event_type{$e_num} eq 'P') {
			print SUM "\t$Otrans{$translate{$event_pos{$e_num}}}\t$Otrans{$translate{$event_pos{$e_num}}}\t@extant\n"; # print Start, End, taxa in the event.
			if ($event_type{$e_num} eq 'S') {
				print ATAB "FT   SNP             $Otrans{$translate{$event_pos{$e_num}}}\nFT                   /colour=1\nFT                   /taxa="; # reading onto original seq
				print TAB "FT   SNP             $rtrans{$translate{$event_pos{$e_num}}}\nFT                   /colour=1\nFT                   /taxa=";
				my $ft = pop(@extant);
				foreach my $t (@extant) {
					print ATAB "taxon_"."$t ";
					print TAB "taxon_"."$t ";
				}
				print TAB "taxon_"."$ft\n";
				print ATAB "taxon_"."$ft\n";
			}
		} elsif ($event_type{$e_num} eq 'R') {
			my @L = split(/\./,$event_pos{$e_num});
			print SUM "\t$Otrans{$translate{$L[0]}}\t$Otrans{$translate{$L[2]}}\t@extant\n";
			foreach my $RecSnp (@{$recsnp{$e_num}}) {
				print SUM "$e_num\tr\t$Otrans{$translate{$RecSnp}}\t$Otrans{$translate{$RecSnp}}\t@extant\n";
			}
			my $RS = scalar(@{$recsnp{$e_num}});
			print ATAB "FT   misc_feature    $Otrans{$translate{$L[0]}}..$Otrans{$translate{$L[2]}}\nFT                   /note=introduces $RS SNPs\nFT                   /taxa=";
			print TAB "FT   misc_feature    $rtrans{$translate{$L[0]}}..$Otrans{$translate{$L[2]}}\nFT                   /note=introduces $RS SNPs\nFT                   /taxa=";
			my $ft = pop(@extant);
			foreach my $t (@extant) {
					print ATAB "taxon_"."$t ";
					print TAB "taxon_"."$t ";
				}
			print TAB "taxon_"."$ft\n";
			print ATAB "taxon_"."$ft\n";
		}
	}
}

#foreach my $position (keys %event_log) {
#	foreach my $type (keys %{$event_log{$position}}) {
#		foreach my $taxon (keys %{$event_log{$position}{$type}}) {
#			if ($type eq "S") {
#				if (defined($event_log{$position}{$type}{$taxon}) && length($event_log{$position}{$type}{$taxon}) > 0 && $event_log{$position}{$type}{$taxon} == 1) {
#					if (defined($rtrans{$position})) {
#						print ATAB "FT   SNP             $Otrans{$translate{$position}}\nFT                   /colour=1\nFT                   /taxa=taxon_"."$taxon\n";
#						print TAB "FT   SNP             $rtrans{$position}\nFT                   /colour=1\nFT                   /taxa=taxon_"."$taxon\n";
#					} else {
#						my $explicable = 0;
#						foreach my $ins_site (keys %log_insertions) {
#							if ($position == 0) {print STDERR "Looking for index $position, position $translate{$position} type $type in taxon $taxon with insertion position $translate{$ins_site}) and extending for $log_insertions{$ins_site}; event binary = $event_log{$position}{$type}{$taxon}\n";}
#							if ($translate{$position} >= $translate{$ins_site-1} && $translate{$position} <= ($translate{$ins_site-1}+$log_insertions{$ins_site})) {
#								print LOG "SNP outside reference in taxon $taxon at $translate{$position} (index $position); inside phage insertion starting at index $ins_site (position $translate{$ins_site-1}) and extending for $log_insertions{$ins_site}\n";
#								$explicable = 1;
#							}
#						}
#						if ($explicable == 0) {
#							print LOG "SNP outside reference in taxon $taxon: index $position, position $translate{$position}; ref site $rtrans{$position}; inexplicable\n";
#						}
#					}
#				}
#			} elsif ($type eq "R") {
#				if (defined($event_log{$position}{$type}{$taxon}) && length($event_log{$position}{$type}{$taxon}) > 0 && $event_log{$position}{$type}{$taxon} == 1) {
#					my @coords = split(/\./,$position);
#					if (defined($rtrans{$coords[0]}) && defined($rtrans{$coords[2]})) {
#						print ATAB "FT   misc_feature    $Otrans{$translate{$coords[0]}}..$Otrans{$translate{$coords[2]}}\nFT                   /taxa=taxon_"."$taxon\n";
#						print TAB "FT   misc_feature    $rtrans{$coords[0]}..$rtrans{$coords[2]}\nFT                   /taxa=taxon_"."$taxon\n";
#						my $Scount = scalar(@{$recsnp{$position}{$type}{$taxon}});
#						if ($#{$recsnp{$position}{$type}{$taxon}} > 1) {
#							print ATAB "FT                   /colour=2\nFT                   /note=introduces $Scount SNPs\n";
#							print TAB "FT                   /colour=2\nFT                   /note=introduces $Scount SNPs\n";
#						} else {
#							print ATAB "FT                   /colour=10\nFT                   /note=introduces $Scount SNP(s)\n";
#							print TAB "FT                   /colour=10\nFT                   /note=introduces $Scount SNP(s)\n";	
#						}
#						foreach my $RecSnp (@{$recsnp{$position}{$type}{$taxon}}) {
#							print TAB "FT   RECSNP          $rtrans{$RecSnp}\nFT                   /colour=4\nFT                   /taxa=taxon_"."$taxon\n";
#							print ATAB "FT   RECSNP          $Otrans{$translate{$RecSnp}}\nFT                   /colour=4\nFT                   /taxa=taxon_"."$taxon\n";
#						}
#					} else {
#						my $explicable = 0;
#						foreach my $ins_site (keys %log_insertions) {
#							if ($translate{$coords[0]} >= $translate{$ins_site-1} && $translate{$coords[2]} <= ($translate{$ins_site-1}+$log_insertions{$ins_site})) {
#								print LOG "Recombination outside reference in taxon $taxon at $translate{$coords[0]} (index $position); inside phage insertion starting at index $ins_site and extending for $log_insertions{$ins_site}\n";
#								$explicable = 1;
#							}
#						}
#						if ($explicable == 0) {
#							print LOG "Recombination outside reference in taxon $taxon: index $coords[0] $coords[2], position $translate{$coords[0]} $translate{$coords[2]}; ref sites $rtrans{$coords[0]} $rtrans{$coords[2]}; inexplicable\n";
#						}
#					}
#				}
#			}
#		}
#	}
#}

close TAB;
close ATAB;

print "I'm done on the simulation \n";

# generate PDF

#system "/nfs/users/nfs_s/sh16/scripts/reportlabtest.py -t $outfile.tree -q taxa -o $outfile.pdf $outfile.tab 2> /dev/null";

# SUBROUTINES

sub croak {
	print STDERR "\ngenerate_taxa.pl -a [alignment file] -m [MGE file] -o [output prefix] -n [# taxa (default: 10)] -b [branching prob (default = 0.05)] -r [rec prob (default = 0.1)] -i [insertion prob (default = 0.1)]\n\n";
}

sub copy_event_log {
	my $t = shift;
	my $i = shift;
	foreach my $en (keys %{$taxon_events{$t}}) {
		if ($taxon_events{$t}{$en} == 1) {
			$taxon_events{$i}{$en} = 1;
			push(@{$event_taxa{$en}},$i);
		}
	}
	
#	foreach my $position (keys %event_log) {
#		foreach my $type (keys %{$event_log{$position}}) {
#			if (defined($event_log{$position}{$type}{$t}) && $event_log{$position}{$type}{$t} == 1) {
#				$event_log{$position}{$type}{$i} = $event_log{$position}{$type}{$t};
#				if ($type eq "R") {
#					@{$recsnp{$position}{$type}{$i}} = @{$recsnp{$position}{$type}{$t}};
#				}
#			}	
#		}
#	}
}

sub mutate {
	my @bases = ('A','C','T','G');
	my $taxon = shift;
	my $seq = $taxa{$taxon};
	my $pm_pos = int(rand(length($seq)-1))+1;
	my $anc_base = uc(substr($seq,$pm_pos,1));
	while ($anc_base eq '-') {							# avoid introducing spurious SNPs into loci where other taxa have phage
		$pm_pos = int(rand(length($seq)-1))+1;
		$anc_base = uc(substr($seq,$pm_pos,1));
	}
	my $der_base = $bases[int(rand(4))];
	while ($der_base eq $anc_base) {
		$der_base = $bases[int(rand(4))];
	}
	substr($seq,$pm_pos,1,$der_base);
	my $Sindex = $back_translate{$pm_pos};
	$exclude{$taxon}{$Sindex} = $Sindex;
	if ($index > 0) {
		$event_log{$e_num} = 1;
		$event_pos{$e_num} = $Sindex;
		$event_type{$e_num} = "S";
		push(@{$event_taxa{$e_num}},$taxon);
		$taxon_events{$taxon}{$e_num} = 1;
		$e_num++;
		#$event_log{$Sindex}{"S"}{$taxon} = 1;
	}
	print LOG "Taxon $taxon mutates to $der_base at position $pm_pos index $Sindex\n";
	return $seq;	
}

sub hom_rec {
	my $taxon = shift;
	my $seq = $taxa{$taxon};
	my @donor_keys = keys %donor;
	my $d_rand = rand @donor_keys; # Choosing a random donor index
	my $donor = $donor_keys[$d_rand]; # Getting that random donor
	my $start = int(rand(length($ancestor)-2))+1;		# phage may have elongated the sequence; must only choose indices from ancestral genome length, not terminal bases to avoid negative recombination lengths
	print STDERR "Start index $start position $translate{$start}:"; # debug
	my $end = $start;
	my $p_end = rand(1);
	while ($p_end > 0.0005) { # Determining how long the recombination event will be,
		$end++;
		$p_end = rand(1);
	}
	if ($end > length($ancestor)) { ## Guarding against super long recombination events!
		$end = length($ancestor) - 1;
	}
	my $excluded = 0; # I think this next bit is working out if the same recombination event has occured previously? Oh no I think its seeing if the proposed recombination is within the bounds of the new sequence (which it should be if aligned to the ancestor)
	foreach my $S (keys %{$exclude{$taxon}}) {
		if (($translate{$start} >= $translate{$S} && $translate{$start} <= $translate{$exclude{$taxon}{$S}}) || ($translate{$S} >= $translate{$start} && $translate{$S} <= $translate{$end}) || ($translate{$start} <= $translate{$S} && $translate{$end} >= $translate{$exclude{$taxon}{$S}}) || ($translate{$start} >= $translate{$S} && $translate{$end} <= $translate{$exclude{$taxon}{$S}})) {
			$excluded = 1;
			print STDERR " - match found to event in $taxon at $translate{$S} to $translate{$exclude{$taxon}{$S}} - ";
		}
	}
	if ($excluded == 0) {
		print STDERR "accepted!\n";
		$exclude{$taxon}{$start} = $end;
		my $length = $end - $start + 1;
		my $donor_seg;
		for (my $c = $start; $c <= $end; $c++) {
			if (defined($log_insertions{$c+1}) && $log_insertions{$c+1} > 0 && $c ne $end) { # This seems odd, not sure what's going on here
				my $dashed = '-' x $log_insertions{$c+1};
				$donor_seg.=$dashed;
			}
			$donor_seg.=substr($donor{$donor},$c,1);
		}
		$length = $translate{$end} - $translate{$start} + 1;
		my $y = length($donor_seg);
		my $old = substr($seq,$translate{$start},$length);
		substr($seq,$translate{$start},$length,$donor_seg);
		if ($index > 0) {
			my $S = $start; my $E = $end;
			my $locus = "$S".".."."$E";
			$event_log{$e_num} = 1;
			$event_pos{$e_num} = $locus;
			$event_type{$e_num} = "R";
			push(@{$event_taxa{$e_num}},$taxon);
			$taxon_events{$taxon}{$e_num} = 1;
			#$event_log{$locus}{"R"}{$taxon} = 1;
			@{$recsnp{$e_num}} = rec_snps($old,$donor_seg,$translate{$start}); # This must be where the recombination snps are calculated
			$e_num++;
			undef($old);
			
			### NEW VERSION ### I currently don't think it's worth keeping track of overwritten mutations
			# as they may be inferred from other sequences in the dataset
			
	#		foreach my $en (keys %{$taxon_events{$taxon}}) {
	#			if ($event_type{$en} eq 'S' && $translate{$event_pos{$en}} >= $S && $translate{$event_pos{$en}} <= $E) {
	#				$taxon_events{$taxon}{$e_num} = 0;			# accounts for recombinations overwriting SNPs
	#			} elsif ($event_type{$en} eq 'R' && ($en+1) ne $e_num) {
	#				my @L = split(/\./,$event_pos{$en});
	#				if ($S <= $L[2] && $E >= $L[0]) {
	#					for (my $idx = 0; $idx <= $#{$recsnp{$en}}; $idx++) {
	#						if ($translate{${$recsnp{$en}}[$idx]} >= $S && $translate{${$recsnp{$en}}[$idx]} <= $E) {
	#							splice(@{$recsnp{$en}},$idx);				# removes recsnps where one recombination overwrites an older one
	#						}
	#					}
	#				}
	#			}
	#		}
		### OLD VERSION ###	
	#		for (my $tpos = $translate{$start};$tpos <= $translate{$end};$tpos++) {		# eliminate the record of SNPs overwritten by recombination
	#			if (defined($event_log{$back_translate{$tpos}}{"S"}{$taxon})) {
	#				$event_log{$back_translate{$tpos}}{"S"}{$taxon} = 0;
	#			}
	#		}
		}
		print LOG "Taxon $taxon gets sequence from $donor between positions $translate{$start} and $translate{$end} (indices: $start and $end)\n";
	} else {
		print STDERR " rejected\n"; # debug	
	}
	return $seq;
}

sub rec_snps {
	my $recip = shift; # The old sequence about to be replaced
	my $donor = shift; # The new sequence to be inserted within
	my $pos = shift; # The start of the recombination event
	my @recip = split(//,$recip); # Create an array from the sequence
	my @donor = split(//,$donor); # Create an array from the donor sequence
	my @snp;

	for (my $i = 0; $i <= $#donor; $i++) { # go from the length of the donor in this for loop
		if ($recip[$i] ne $donor[$i] && $recip[$i] ne '-' && $donor[$i] ne '-') { # This takes into account the gaps present in donor and recipient and not counting them as SNPs, so where is this going wrong
				push(@snp,$back_translate{($i+$pos+1)}); # Is the original sequence just really divergent to the one we see in the alignment?
		}
	}

	return @snp;
}

sub phage_ins {
	my $taxon = shift;
	my $seq = $taxa{$taxon};
	my @mge_keys = keys %mge;
	my $m_rand = rand @mge_keys;
	my $mge_name = $mge_keys[$m_rand];
	my $ins_site = int(rand(length($ancestor)-3))+2;		# ensure insertion site exists in the ancestor, and for practical indexing reasons make sure these events do not occur near the initial or terminal bases
	while (defined($superinfection{$taxon}{$ins_site}) && $superinfection{$taxon}{$ins_site} == 1) {
		$ins_site = int(rand(length($ancestor)-3))+2;
	}
	$superinfection{$taxon}{$ins_site} = 1;
	my $ins_position = $translate{$ins_site};
	my $dashed = '-' x length($mge{$mge_name});
	
	# add phage to sequences
	foreach my $t (keys %taxa) {
		my $before = substr($taxa{$t},0,($translate{$ins_site}-1));
		my $end_length = length($taxa{$t}) - $translate{$ins_site};
		my $after = substr($taxa{$t},($translate{$ins_site}-1));
		my $l = length($taxa{$t});
		print LOG "Taxon $t ($taxon) length $l pos $translate{$ins_site}\n";
		if ($t eq $taxon) {
			$taxa{$t} = $before.$mge{$mge_name}.$after;
		} else {
			$taxa{$t} = $before.$dashed.$after;
		}
	}
	
	# increment translation table following phage insertion, new indices for translation
	foreach my $i (sort keys %translate) {
		if ($translate{$i} > $translate{$ins_site}) {
			$translate{$i}+=length($mge{$mge_name});
		}
	}
	$translate{$ins_site}+=length($mge{$mge_name});
	$i_max++;
	for (my $n = $translate{$ins_site-1}+1; $n < ($ins_position+length($mge{$mge_name}));$n++) {
		$translate{$i_max} = $n;
		$i_max++;
	}
	my $t_ins = $back_translate{$ins_site};
	$log_insertions{$ins_site}+=length($mge{$mge_name});
	
	$event_log{$e_num} = 1;
	$event_type{$e_num} = 'P';
	$event_pos{$e_num} = $ins_site;
	push(@{$event_taxa{$e_num}},$taxon);
	$taxon_events{$taxon}{$e_num} = 1;
	$e_num++;
	
	back_translate();
	print LOG "Taxon $taxon acquires phage $mge_name at position $ins_site (insertion length $log_insertions{$ins_site}) indexed at $t_ins\n";	
}

sub back_translate {
	map {$back_translate{$translate{$_}} = $_;} (keys %translate);
}
