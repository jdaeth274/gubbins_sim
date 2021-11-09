#! /usr/bin/perl
use warnings;
use Bio::AlignIO;
use Getopt::Long;

# parse input options

my $aln_file;
my $outfile;
my $true_file;
my $gub_file;
my $gub_snp_file;
my $help = 0;

GetOptions (
	"h" => \$help,
	"a=s" => \$aln_file,
	"o=s" => \$outfile,
	"t=s" => \$true_file,
	"g=s" => \$gub_file,
	"s=s" => \$gub_snp_file
);

if ($help == 1) {
	croak(); exit(0);
}


unless (defined($aln_file) && defined($outfile) && defined($true_file) && defined($gub_file)) {
	print STDERR "Undefined file!\n";print STDERR "Error: $!\n" and croak(); exit(1);
}

# parse alignment, record names of taxa and length of alignment

my $alnObject = Bio::AlignIO->new(-file => $aln_file, -format => "Fasta");
my $aln = $alnObject->next_aln;

my @taxa;

foreach my $seq ($aln->each_seq) {
	my $id = $seq->id;
	push(@taxa,$id);
}

# find 'true' recombinations

open TRUE, $true_file or print STDERR "Error: $!\n" and croak() and exit(1);

my %snppos;
my %snptaxa;
my %recsnppos;
my %recsnptaxa;
my %recstart;
my %recend;
my %rectaxa;
my $rindex = 0;
my %branchsnps;
my %inrecsnps;
my $total_true_snps = 0;

foreach (<TRUE>) {
	unless (/^Event/) {
		chomp;
		my @F = split(/\s+/,$_);
		if ($F[1] eq 'S') {
			$snppos{$F[0]} = $F[2];
			my @d = split(/,/,$F[4]);
			@{$snptaxa{$F[0]}} = @d;
			$total_true_snps++;
		} elsif ($F[1] eq 'R') {
			$recstart{$F[0]} = $F[2];
			$recend{$F[0]} = $F[3];
			my @d = split(/,/,$F[4]);
			@{$rectaxa{$F[0]}} = @d;
		} elsif ($F[1] eq 'r') {
			$recsnppos{$rindex} = $F[2];
			my @d = split(/,/,$F[4]);
			@{$recsnptaxa{$rindex}} = @d;
			push(@{$inrecsnps{$F[0]}},$rindex);
			$total_true_snps++;
			$rindex++;
		}
	}
}

my $num_branches = scalar(@taxa)*2-3;
my $snps_per_branch = $total_true_snps/$num_branches;

close TRUE;

# find predicted recombinations

open GUB, $gub_file or print STDERR "Error: $!\n" and croak() and exit(1);

$in_rec = 0;
my %pred_rec;
my %prec_start;
my %prec_end;
my $prec = 0;

foreach (<GUB>) {
	chomp;
	if (/^FT   misc_feature/) {
		my @F = split(/\s+/,$_);
		$locus = $F[2];
		$in_rec = 1;
	} elsif ($in_rec == 1 && /FT                   \/taxa/) {
		$_ =~ s/\"//g;
		$_ =~ s/taxon_//g;
		my @data = split(/=/,$_);
		my @F = split(/, /,$data[1]);
		@{$pred_rec{$prec}} = @F;
		my @coords = split(/\./,$locus);
		$prec_start{$prec} = $coords[0];
		$prec_end{$prec} = $coords[2];
		$prec++;
	}
}

close GUB;

# filter predicted recombinations to get rid of duplicates

foreach my $locA (keys %prec_start) {
	foreach my $locB (keys %prec_start) {
		unless ($locA eq $locB) {
			if (compare_arrays($pred_rec{$locA},$pred_rec{$locB})) {
				if ($prec_end{$locB} <= $prec_end{$locA} && $prec_start{$locB} >= $prec_start{$locA}) {
					delete($pred_rec{$locB});
					delete($prec_start{$locB});
					delete($prec_end{$locB});
				}
			}
		}
	}
}

# find SNP reconstruction

open SNP, $gub_snp_file or print STDERR "Error: $!" and croak() and exit(1);

my %p_snppos;
my %p_snptaxa;
my $snp_i = 0;
my $pos;
my $pred_rec_snp = 0;

foreach (<SNP>) {
	chomp;
	if (/^FT   SNP/) {
		my @F = split(/\s+/,$_);
		$pos = $F[2];
	} elsif (/^FT                   \/taxa/) {
		$_ =~ s/\"//g;
		$_ =~ s/taxon_//g;
		my @data = split(/=/,$_);
		my @F = split(/\s+/,$data[1]);
		$p_snppos{$snp_i} = $pos;
		@{$p_snptaxa{$snp_i}} = @F;
		$snp_i++;
	}
}

close SNP;

# which snps are in which predicted recombinations

my %p_inrecsnps;

foreach $snp_i (keys %p_snppos) {
	foreach $prec (keys %pred_rec) {
		my @F = ($prec_start{$prec},$prec_end{$prec});
		if ($p_snppos{$snp_i} >= $prec_start{$prec} && $p_snppos{$snp_i} <= $prec_end{$prec}) {
			my @compA = sort(@{$pred_rec{$prec}}); my @compB = sort(@{$p_snptaxa{$snp_i}});
			if (compare_arrays(\@compA,\@compB)) {
				push(@{$p_inrecsnps{$prec}},$snp_i);
				$pred_rec_snp++;
			}
		}
	}
}

# find stats of truth and prediction

my $true_rec = 0;
my $true_length = 0;
my $pred_rec = 0;
my $pred_length;

foreach $locus (keys %recstart) {
	$true_rec++;
	$true_length+=(1+$recend{$locus}-$recstart{$locus});
}

foreach $locus (keys %pred_rec) {
	$pred_rec++;
	$pred_length+=(1+$prec_end{$locus} - $prec_start{$locus});
}

# compare predictions with the truth

my $FN_length = 0;
my $FP_length = 0;
my $TP_length = 0;
my $FN_event = 0;
my $FP_event = 0;
my $TP_event = 0;
my $FN_snp = 0;
my $FP_snp = 0;
my $TP_snp = 0;
my $in_FP_event = 0;

my %found_event; map {$found_event{$_} = 0;} (keys %pred_rec);
my $Lbound; my $Rbound;

# first, establish links between predicted and true recombinations (may be one to one or many to many)

my $link = 0;
my $max_link = 0;
my %truelink;
my %plinkback;
my %tlinkback;
my %predlink;

foreach my $locT (keys %recstart) {
	foreach my $locP (keys %pred_rec) {
		if (($recstart{$locT} >= $prec_start{$locP} && $recstart{$locT} <= $prec_end{$locP}) || ($prec_start{$locP} >= $recstart{$locT} && $prec_start{$locP} <= $recend{$locT})) {
			if (check_taxa($rectaxa{$locT},$pred_rec{$locP})) {
				my $linkA; my $linkB;
				if (defined($plinkback{$locP})) {
					$linkA = $plinkback{$locP};
				}
				if (defined($tlinkback{$locT})) {
					$linkB = $tlinkback{$locT};
				}
				if (!(defined($linkA)) && !defined($linkB)) {
					$link = $max_link = $max_link+1;
				} elsif (defined($linkA) && defined($linkB) && $linkA ne $linkB) {
					foreach my $locX (keys %{$truelink{$linkB}}) {
						$truelink{$linkA}{$locX} = 1;
						$tlinkback{$locX} = $linkA;
					}
					undef(%{$truelink{$linkB}});
					foreach my $locX (keys %{$predlink{$linkB}}) {
						$predlink{$linkA}{$locX} = 1;
						$plinkback{$locX} = $linkA;
					}
					undef(%{$predlink{$linkB}});
					$link = $linkA;
				} elsif (defined($linkA)) {
					$link = $linkA;
				} elsif (defined($linkB)) {
					$link = $linkB;
				}
				$truelink{$link}{$locT} = 1;
				$predlink{$link}{$locP} = 1;
				$plinkback{$locP} = $link;
				$tlinkback{$locT} = $link;
			}
		}
	}
}

# all true recombinations without links are false negatives

foreach my $locT (keys %recstart) {
	unless (defined($tlinkback{$locT}) && scalar(keys %{$predlink{$tlinkback{$locT}}}) > 0) {
		$FN_event++;
		$FN_length+=(1+$recend{$locT}-$recstart{$locT});
		if (defined(@{$inrecsnps{$locT}})) {
			$FN_snp+=scalar(@{$inrecsnps{$locT}});
		}
	}
}

# all predicted recombinations without links are false positives

foreach my $locP (keys %pred_rec) {
	unless (defined($plinkback{$locP}) && scalar(keys %{$truelink{$plinkback{$locP}}}) > 0) {
		$FP_event++;
		$FP_length+=(1+$prec_end{$locP}-$prec_start{$locP});
		$FP_snp+=scalar(@{$p_inrecsnps{$locP}});
		$in_FP_event+=scalar(@{$p_inrecsnps{$locP}});
	}
}

# where links exist, work out the level of overlap

foreach $link (keys %truelink) {
	if (scalar(keys %{$truelink{$link}}) > 0 && scalar(keys %{$predlink{$link}}) > 0) {			# this gets rid of links that were deleted when merging clusters of recombinations
		my $snpP; my $snpT; my $Plen; my $Tlen;
		foreach my $locT (keys %{$truelink{$link}}) {
			$Tlen+=(1+$recend{$locT}-$recstart{$locT});
			if (defined(@{$inrecsnps{$locT}})) {
				push(@{$snpT},@{$inrecsnps{$locT}});
			}
			$TP_event++;
		}
		foreach my $locP (keys %{$predlink{$link}}) {
			$Plen+=(1+$prec_end{$locP}-$prec_start{$locP});
			push(@{$snpP},@{$p_inrecsnps{$locP}});
		}
		my $tp = 0;
		foreach my $locT (keys %{$truelink{$link}}) {
			my @F = ($recstart{$locT},$recend{$locT});
			foreach my $locP (keys %{$predlink{$link}}) {
				my @G = ($prec_start{$locP},$prec_end{$locP});
				if ($G[0] >= $F[0] && $G[1] <= $F[1]) {	# prediction is within real
					$tp+=(1+$G[1]-$G[0]);	
				} elsif ($F[0] >= $G[0] && $F[1] <= $G[1]) {	# real is within predicted
					$tp+=(1+$F[1]-$F[0]);
				} elsif ($G[0] <= $F[0] && $G[1] <= $F[1] && $G[1] > $F[0]) {	# overlap at start
					$tp+=(1+$G[1]-$F[0]);
				} elsif ($G[0] >= $F[0] && $F[1] <= $G[1] && $G[0] <= $F[1]) {	# overlap at end
					$tp+=(1+$F[1]-$G[0]);
				}
			}
		}
		#unless (defined(@{$snpT})) {print "not defined: @{$snpT} link $link\n";foreach my $locT (keys %{$truelink{$link}}) {print STDERR "locT $locT s $recstart{$locT} e $recend{$locT} SNPs @{$inrecsnps{$locT}}\n";} foreach my $locP (keys %{$predlink{$link}}) {print STDERR "locP $locP s $prec_start{$locP} e $prec_end{$locP} SNPs @{$p_inrecsnps{$locP}}\n";} die;}
		if (defined(@{$snpT}) && $#{$snpT} > -1) {
			new_eval_rec_snps($snpT,$snpP);
		} else {
			$FP_snp+=scalar(@{$snpP});	# this problem arises when a zero SNP recombination occurs on the earliest branch and it's reconstructed as introduction of reference sequence by a recombination
		}
		$TP_length+=$tp;
		$FN_length+=($Tlen-$tp);
		$FP_length+=($Plen-$tp);
	}
}

# print output

open OUT, "> $outfile" or print STDERR "Error: $!\n" and croak() and exit(1);

print OUT "Alignment\t#True\t#Pred\tTP_events\tFP_events\tFN_events\tLen(True)\tLen(Pred)\tTP_length\tFP_length\tFN_length\t#TrueSNP\t#PredSNP\tSNPperBranch\tTP_snp\tFP_snp\tFN_snp\tInFPEvent\n";
print OUT "$aln_file\t$true_rec\t$pred_rec\t$TP_event\t$FP_event\t$FN_event\t$true_length\t$pred_length\t$TP_length\t$FP_length\t$FN_length\t$rindex\t$pred_rec_snp\t$snps_per_branch\t$TP_snp\t$FP_snp\t$FN_snp\t$in_FP_event\n";

# SUBROUTINES

sub new_eval_rec_snps {
	my $tp = 0; my $fp = 0; my $fn = 0;
	### here need snp reconstruction!
	my @Tr = @{$_[0]};
	my @Pr = @{$_[1]};
	foreach my $T (@Tr) {
		my $match = 0; my $P = 0;
		while ($P <= $#Pr) {
			if ($p_snppos{$Pr[$P]} == $recsnppos{$T}) {
				$match = 1;
				splice(@Pr,$P,1);
				last;
			} else {
				$P++;
			}
		}
		if ($match == 1) {
			$tp++;
		} else {
			$fn++;
		}
	}
	$fp = scalar(@Pr);
	$FP_snp+=$fp;
	$TP_snp+=$tp;
	$FN_snp+=$fn;
}

sub check_taxa {
	my @actual = sort(@{$_[0]});
	my @pred = sort(@{$_[1]});
	my %p; map {$p{$_} = 1;} @pred;
	my @unpred;
	foreach my $t (@taxa) {
		$t =~s/taxon_//g;
		unless (defined($p{$t}) && $p{$t} == 1) {
			push(@unpred,$t);
		}
	}
	my @u = sort(@unpred);
	if (compare_arrays(\@actual,\@pred)) {
		return 1; 
	} elsif (compare_arrays(\@actual,\@u)) {
		return 1; 
	} else {
		return 0;
	}
	
}

sub compare_arrays {
	my ($first, $second) = @_;
	no warnings; # silence spurious -w undef complaints
	return 0 unless @{$first} == @{$second};
	for (my $i = 0; $i < @{$first}; $i++) {
		return 0 if $first->[$i] ne $second->[$i];
	}
	return 1;
}

sub croak {
	print STDERR "Usage: gubbins_evaluation.pl -t [simulation summary file] -a [input alignment] -g [gubbins recombination tab] -s [gubbins SNP tab] -o [outfile]\n";
}
