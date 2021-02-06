#!/usr/bin/perl
use strict; use warnings;
use FAlite;
use DataBrowser;
use File::Basename;
use Getopt::Std;
our ($opt_t, $opt_i, $opt_I, $opt_e, $opt_E, $opt_m, $opt_x, $opt_s);
getopts('ti:I:e:E:m:x:s');

my $MIN_INTRON = 30;
my $MAX_INTRON = 1000;
my $MIN_EXON   = 15;
my $MAX_EXON   = 1000;
my $ORDER      = 0;
my $EXPLICIT   = 0;

die "
usage: $0 <ann> <dna> <don> <acc> <gc> <inter> <title>
  <ann>   annotation file in ZFF
  <dna>   sequence file in FASTA
  <don>   length of intronic portion of donor site
  <acc>   length of intronic portion of acceptor site
  <gc>    GC content of intergenic portion (eg. 0.4)
  <inter> intergenic fraction of genome (eg. 0.7)
  <title> something unique that includes the genome name
options:
  -t         test mode, keep temp directory        [default remove]
  -i <int>   minimum intron length                 [$MIN_INTRON]
  -I <int>   maximum intron length                 [$MAX_INTRON]
  -e <int>   minimum exon length                   [$MIN_EXON]
  -E <int>   maximum exon length                   [$MAX_EXON]
  -m <int>   use Markov models of <int> order      [$ORDER]
  -x <int>   use explict durations of <int> length [default off]
  -s         keep single exon genes                [default skip]
" unless @ARGV == 7;
my ($ANN, $DNA, $DON, $ACC, $GC, $INTER, $TITLE) = @ARGV;
die "GC content out of range\n" if $GC < 0.2 or $GC > 0.8;
die "Intergenic fraction out of range\n" if $INTER < 0.01 or $INTER > 0.99;
my $AT = 1 - $GC;

$MIN_INTRON = $opt_i if $opt_i;
$MAX_INTRON = $opt_I if $opt_I;
$MIN_EXON   = $opt_e if $opt_e;
$MAX_EXON   = $opt_E if $opt_E;
$ORDER      = $opt_m if $opt_m;
$EXPLICIT   = $opt_x if $opt_x;
my $KEEPTEMP = $opt_t ? $opt_t : 0;
my $KEEPSNGL = $opt_s ? $opt_s : 0;

my $MER = $ORDER +1;

my %E; # exon
my %I; # intron
my @D; # donor
my @A; # acceptor
my @LE; # length exon
my @LI; # length intron

# setup
my $DIR = "temp.training.$$";
END {
	if (defined $DIR and not $KEEPTEMP) {
		chdir "..";
		system("rm -rf $DIR");
	}
}
mkdir $DIR;
my $ann_file = basename($ANN);
my $dna_file = basename($DNA);
run("cp $ANN $DIR/$ann_file");
run("cp $DNA $DIR/$dna_file");
chdir $DIR;

# skip single exon genes (because they don't splice)
if ($KEEPSNGL) {
	run("ln -s $ann_file select.ann");
	run("ln -s $dna_file select.dna");
} else {
	my @skip;
	if ($ann_file =~ /gz$/) {open(IN, "gunzip -c $ann_file |") or die}
	else                    {open(IN, $ann_file) or die}
	my $id;
	while (<IN>) {
		if (/^>(\S+)/) {$id = $1}
		elsif (/^Esngl/) {push @skip, $id}
	}
	close IN;
	open(OUT, ">list") or die;
	print OUT join("\n", @skip), "\n";
	close OUT;
	run("data-selector.pl -v $ann_file $dna_file list");
}

# require plus strand genes in the uni, alt, and wrn classes, skips olp, err
run("fathom -categorize 3 select.ann select.dna");
run("cat uni.ann alt.ann wrn.ann > genes.ann");
run("cat uni.dna alt.dna wrn.dna > genes.dna");
run("fathom -export 3 -plus genes.ann genes.dna");
run("fathom -exon-intron export.ann export.dna > exon-intron");

# training loop
open(IN, "exon-intron") or die;
while (<IN>) {
	chomp;	
	if    (/^[a-z]/) {trainIntron(uc $_)}
	elsif (/^[A-Z]/) {trainExon(uc $_)}
}
close IN;

# make sure there are no zero counts in the exon or introns
my $slots = 4 ** $ORDER;
my $e_ctx = keys %E;
my $i_ctx = keys %I;
die "Exon model order too high\n"   if $ORDER and $e_ctx != $slots;
die "Intron model order too high\n" if $ORDER and $i_ctx != $slots;

# convert counts to probabilities
count_to_prob(\%E);
count_to_prob(\%I);
count_to_prob2(\@A);
count_to_prob2(\@D);

my $Estats = seqstats(\@LE);
my $Istats = seqstats(\@LI);

my $ExonLength   = $Estats->{average} -3 -3;
my $IntronLength = $Istats->{average} -$DON -$ACC;

my $Nfrac = $INTER;
my $genic = 1 - $INTER;
my $Total = $Estats->{total} + $Istats->{total};
my $Efrac = $genic * $Estats->{total} / $Total;
my $Ifrac = $genic * $Istats->{total} / $Total;

#############
# write HMM #
#############

my ($initN, $termN) = ($Nfrac, $Nfrac);
my ($initI, $termI) = ($Ifrac, $Ifrac);
my ($initE, $termE) = ($Efrac/2, $Efrac/2);

my $leaveE = 1/$ExonLength;
my $stayE  = 1 - $leaveE;
my $leaveI = 1/$IntronLength;
my $stayI  = 1 - $leaveI;
my $STATES = 1 + 1 + 3+$DON + 1 + $ACC+3 + 1; # N,E,D,I,A,E

# header
print "# denada HMM for $TITLE splicing\n";
print "# $STATES states\n"; # explicit mention when it happens
printf "# Intergenic: %.3f GC, %.3f of genome\n", $GC, $Nfrac;
printf "# Exon: %.3f bp average\n", $ExonLength;
printf "# Intron: %.3f bp average\n", $IntronLength;
print "# Markov: order $ORDER (applies to Exon & Intron only)\n";
print "# Donor: 3 nt exon + $DON nt intron\n";
print "# Acceptor: $ACC nt intron + 3 nt exon\n";
print "# Filter: $MIN_EXON < exon > $MAX_EXON, $MIN_INTRON < intron > $MAX_INTRON\n\n";

print "denada-HMM $TITLE-splicing $STATES\n\n";

# intergenic state
printf "denada-State Inter %.3f %.3f 1 0 0\n", $initN, $termN;
print "\tInter 1.0\n";
printf "\t%.3f\t%.3f\t%.3f\t%.3f\n\n", $AT/2, $GC/2, $GC/2, $AT/2;

# exon 1
printf "denada-State Exon1 %.3f %.3f 2 $ORDER 0\n", $initE, $termE;
printf "\tExon1 %.3f\n\tD1    %.3f\n", $stayE, $leaveE;
foreach my $ctx (sort keys %E) {
	printf "\t%.3f\t%.3f\t%.3f\t%.3f\n", $E{$ctx}{A}, $E{$ctx}{C}, $E{$ctx}{G}, $E{$ctx}{T};
}
print "\n";

# donor
for (my $i = 0; $i < @D -1; $i++) {
	printf "denada-State D%d 0.0 0.0 1 0 0\n\tD%d 1.0\n", $i+1, $i+2;
	printf "\t%.3f\t%.3f\t%.3f\t%.3f\n\n", $D[$i]{A}, $D[$i]{C}, $D[$i]{G}, $D[$i]{T};
}
printf "denada-State D%d 0.0 0.0 1 0 0\n\tIntron 1.0\n", scalar(@D);
printf "\t%.3f\t%.3f\t%.3f\t%.3f\n\n", $D[@D-1]{A}, $D[@D-1]{C}, $D[@D-1]{G}, $D[@D-1]{T};

# intron
printf "denada-State Intron %.3f %.3f 2 $ORDER 0\n", $initI, $termI;
printf "\tIntron %.3f\n\tA1     %.3f\n", $stayI, $leaveI;
foreach my $ctx (sort keys %I) {
	printf "\t%.3f\t%.3f\t%.3f\t%.3f\n", $I{$ctx}{A}, $I{$ctx}{C}, $I{$ctx}{G}, $I{$ctx}{T};
}
print "\n";

# acceptor
for (my $i = 0; $i < @A -1; $i++) {
	printf "denada-State A%d 0.0 0.0 1 0 0\n\tA%d 1.0\n", $i+1, $i+2;
	printf "\t%.3f\t%.3f\t%.3f\t%.3f\n\n", $A[$i]{A}, $A[$i]{C}, $A[$i]{G}, $A[$i]{T};
}
printf "denada-State A%d 0.0 0.0 1 0 0\n\tExon2 1.0\n", scalar(@A);
printf "\t%.3f\t%.3f\t%.3f\t%.3f\n\n", $A[@A-1]{A}, $A[@A-1]{C}, $A[@A-1]{G}, $A[@A-1]{T};

# exon 2
printf "denada-State Exon2 %.3f %.3f 1 $ORDER 0\n", $initE, $termE;
printf "\tExon1 %.3f\n", $stayE;
foreach my $ctx (sort keys %E) {
	printf "\t%.3f\t%.3f\t%.3f\t%.3f\n", $E{$ctx}{A}, $E{$ctx}{C}, $E{$ctx}{G}, $E{$ctx}{T};
}
print "\n";


###############
# subroutines #
###############

sub run {
	my ($cmd) = @_;
	system($cmd) == 0 or die "ERROR: running $cmd\nTo check intermediate files use -t\n";
}

sub trainExon {
	my ($seq) = @_;
		
	# filter & length
	my $length = length($seq);
	return if ($length < $MIN_EXON or $length > $MAX_EXON);
	push @LE, $length;
	
	# sequence
	for (my $i = 3; $i < length($seq) - $MER -3; $i++) {
		my $ctx = substr($seq, $i, $MER -1);
		my $nt  = substr($seq, $i + $MER -1, 1);
		next if $ORDER and $ctx !~ /^[ACGT]+$/;
		next if $nt !~ /^[ACGT]$/;
		$E{$ctx}{$nt}++;
	}
	
	# splice donor
	for (my $i = 0; $i < 3; $i++) {
		my $nt = substr($seq, $i, 1);
		next if $nt !~ /^[ACGT]$/;
		$D[$i]{$nt}++;
	}
	
	# splice acceptor
	for (my $i = 0; $i < 3; $i++) {
		my $nt = substr($seq, $i, 1);
		next if $nt !~ /^[ACGT]$/;
		$A[$i+$ACC]{$nt}++;
	}
	
}

sub trainIntron {
	my ($seq) = @_;
	
	# filter & length
	my $length = length($seq);
	return if ($length < $MIN_INTRON or $length > $MAX_INTRON);
	push @LI, $length;
	
	# sequence
	for (my $i = 3; $i < length($seq) - $MER -3; $i++) {
		my $ctx = substr($seq, $i, $MER -1);
		my $nt  = substr($seq, $i + $MER -1, 1);
		next if $ORDER and $ctx !~ /^[ACGT]+$/;
		next if $nt !~ /^[ACGT]$/;
		$I{$ctx}{$nt}++;
	}
	
	# splice donor
	for (my $i = 0; $i < $DON; $i++) {
		my $nt = substr($seq, $i, 1);
		next if $nt !~ /^[ACGT]$/;
		$D[$i+3]{$nt}++;
	}
	
	# splice acceptor
	my $acc = substr($seq, -$ACC);
	for (my $i = 0; $i < $ACC; $i++) {
		my $nt = substr($acc, $i, 1);
		next if $nt !~ /^[ACGT]$/;
		$A[$i]{$nt}++;
	}
}


sub count_to_prob2 {
	my ($data) = @_;
	my @nt = qw(A C G T);
	for (my $i = 0; $i < @$data; $i++) {
		my $total = 0;
		foreach my $nt (@nt) {$data->[$i]{$nt} = 0 if not defined $data->[$i]{$nt}}
		foreach my $nt (@nt) {$total += $data->[$i]{$nt}}
		foreach my $nt (@nt) {$data->[$i]{$nt} /= $total}
	}
}

sub count_to_prob {
	my ($data) = @_;
	foreach my $ctx (keys %$data) {
		my $total = 0;
		foreach my $nt (keys %{$data->{$ctx}}) {$total += $data->{$ctx}{$nt}}
		foreach my $nt (keys %{$data->{$ctx}}) {$data->{$ctx}{$nt} /= $total}
	}
}


sub seqstats {
	my ($ary) = @_;
	
	my @len = sort {$a <=> $b} @$ary;
	my $mode = $len[@len/2]; # considering mode rather than average...
	
	my $total = 0;
	foreach my $l (@len) {$total += $l}
	my $ave = $total/@len;
	
	return {
		average => $ave,
		mode => $mode,
		total => $total
	};
}
