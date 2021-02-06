#!/usr/bin/perl
use strict; use warnings;
use IME;
use FAlite;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_w $opt_m $opt_s $opt_p);
getopts('w:s:p:');

my $WORDSIZE  = 5;
my $IME_START = 100;
my $PROXIMAL  = 500;
my $MIN_GENE  = 1000;

die "
usage: imeter-hmm-trainer.pl <genome file> <exon-intron file>
options:
  -w <int>  wordsize [$WORDSIZE]
  -s <int>  IME start [$IME_START]
  -m <int>  minimum gene length [$MIN_GENE]
  -p <int>  proximal distance [$PROXIMAL]
" unless @ARGV == 2;
my ($GENOME, $EXON_INTRON) = @ARGV;

$WORDSIZE  = $opt_w if $opt_w;
$IME_START = $opt_s if $opt_s;
$PROXIMAL  = $opt_p if $opt_p;
$MIN_GENE  = $opt_m if $opt_m;

my ($Exon, $IME, $Intron, $Genome) = ({}, {}, {}, {});
open(EI, $EXON_INTRON) or die;
my $fasta = new FAlite(\*EI);
while (my $entry = $fasta->nextEntry) {
	my $seq = $entry->seq;
	next unless length($seq) > $MIN_GENE;
	my $dist = substr($seq, $PROXIMAL);
	next unless defined $dist and length($dist) >= $WORDSIZE;
	
	# count exons
	my @exon = split(/[a-z]+/, $dist);
	foreach my $exon (@exon) {count($Exon, $exon)}
	
	# count introns
	my @intron= split(/[A-Z]+/, $dist);
	foreach my $intron (@intron) {count($Intron, uc $intron)}
	
	# count proximal region
	my $prox = uc substr($seq, $IME_START, $PROXIMAL - $IME_START);
	count($IME, $prox);
	
}
close EI;

open(IN, "gunzip -c $GENOME |") or die;
$fasta = new FAlite(\*IN);
while (my $chrom = $fasta->nextEntry) {
	count($Genome, $chrom->seq); # I suppose I could reverse-complement...
}
close IN;

freq($Exon);
freq($Intron);
freq($IME);
freq($Genome);

# 230, 157

my $ctx = $WORDSIZE - 1;
print "denada-HMM proximeter-w$WORDSIZE-s$IME_START-p$PROXIMAL-m$MIN_GENE 4

denada-State Gen 0.75 0.75 4 $ctx 0\n
	Gen  0.997
	Exon 0.001
	Int  0.001
	Prox 0.001\n
", table($Genome), "
denada-State Exon 0.10 0.10 4 $ctx 0\n
	Gen  0.002
	Exon 0.995
	Int  0.002 
	Prox 0.001\n
", table ($Exon), "
denada-State Int 0.10 0.10 4 $ctx 0\n
	Gen  0.001
	Exon 0.003
	Int  0.995
	Prox 0.001\n
", table ($Intron), "
denada-State Prox 0.05 0.05 4 $ctx 0\n
	Gen  0.001
	Exon 0.001
	Int  0.001
	Prox 0.997\n
", table ($IME), "\n";

###############################################################################
# subroutines
###############################################################################

sub table {
	my ($data) = @_;
	my $out = "";
	foreach my $ctx (sort keys %$data) {
		$out .= "\t";
		foreach my $nt ('A', 'C', 'G', 'T') {
			my $val = sprintf "%.4f ", $data->{$ctx}{$nt};
			$out .= $val
		}
		$out .= "\n";
	}
	return $out;
}

sub freq {
	my ($data) = @_;
	foreach my $ctx (keys %$data) {
		my $total = 0;
		foreach my $nt ('A', 'C', 'G', 'T') {
			$total += $data->{$ctx}{$nt};
		}
		foreach my $nt ('A', 'C', 'G', 'T') {
			$data->{$ctx}{$nt} /= $total;
		}
	}
}

sub count {
	my ($data, $seq) = @_;
	return unless length($seq) >= $WORDSIZE;
	for (my $i = 0; $i < length($seq) - $WORDSIZE + 1; $i++) {
		my $word = substr($seq, $i, $WORDSIZE);
		next unless $word =~ /^[ACGT]+$/i;
		my $ctx = substr($word, 0, $WORDSIZE -1);
		my $nt = substr($word, -1, 1);
		$data->{$ctx}{$nt}++;
	}
}

__END__

still uses the dreaded exon-intron file...
more of a proof-of-concept than something to be used
