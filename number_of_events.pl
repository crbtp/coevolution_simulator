#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use 5.010;
use List::Util qw(min max);
use Getopt::Long;

#perl do.pl -r 0.1 -f 0.5
my ($rate,$fitp);
GetOptions(
        'r|rate=f'  => \$rate,
        'f|fitp=f'  => \$fitp,
);


my $child=10;#number of offspring
my $generation=1000;#number of generation
my $individual=1000;#expected number of individual


my %hashA;
my %hashB;
my %hashC;

#A is origin AA, C is coevolved AA, B is others
$hashA{"A"}=1-$rate;
$hashA{"B"}=1-$rate+0.95*$rate;
$hashA{"C"}=1;

$hashB{"B"}=(1-$rate)+0.9*$rate;
$hashB{"A"}=(1-$rate)+0.9*$rate+0.05*$rate;
$hashB{"C"}=1;

$hashC{"C"}=1-$rate;
$hashC{"B"}=1-$rate+0.95*$rate;
$hashC{"A"}=1;

my @seq=qw(AA);
my @nextseq=qw();
my $n=1;
my $mark=0;#total mutation
my $mark1=0;#successive mutation
my $mark2=0;#double mutation
my @probA;
my @probC;
while ($n<=$generation) {
	##########################genarate offspring#############################
	foreach my $key (@seq) {
		my $n1=1;
		while ($n1<=$child) {
			$n1++;
			$key=~m/(.)(.)/;
			my $p1=$1;
			my $p2=$2;

			################mutation 1###################
			my $j=rand(1);
			my $np1;
			if ($p1 eq "A") {
				if ($j<=$hashA{"A"}) {
					$np1="A";
				}
				elsif ($j<=$hashA{"B"}) {
					$np1="B";
				}
				else {
					$np1="C";
				}
			}
			if ($p1 eq "B") {
				if ($j<=$hashB{"B"}) {
					$np1="B";
				}
				elsif ($j<=$hashB{"A"}) {
					$np1="A";
				}
				else {
					$np1="C";
				}
			}
			if ($p1 eq "C") {
				if ($j<=$hashC{"C"}) {
					$np1="C";
				}
				elsif ($j<=$hashC{"B"}) {
					$np1="B";
				}
				else {
					$np1="A";
				}
			}
			################mutation 2###################
			$j=rand(1);
			my $np2;
			if ($p2 eq "A") {
				if ($j<=$hashA{"A"}) {
					$np2="A";
				}
				elsif ($j<=$hashA{"B"}) {
					$np2="B";
				}
				else {
					$np2="C";
				}
			}
			if ($p2 eq "B") {
				if ($j<=$hashB{"B"}) {
					$np2="B";
				}
				elsif ($j<=$hashB{"A"}) {
					$np2="A";
				}
				else {
					$np2="C";
				}
			}
			if ($p2 eq "C") {
				if ($j<=$hashC{"C"}) {
					$np2="C";
				}
				elsif ($j<=$hashC{"B"}) {
					$np2="B";
				}
				else {
					$np2="A";
				}
			}
			############################################
			my $new=$np1.$np2;
			push @nextseq,$new;
			if ($p1 ne $np1) {
				$mark++;
			}
			if ($p2 ne $np2) {
				$mark++;
			}
			if ($key eq "AC" or $key eq "CA" or $key eq "BC" or $key eq "CB") {
				if ($new eq "CC") {
					$mark1++;
				}
			}
			if ($key eq "AB" or $key eq "BA" or $key eq "AA" or $key eq "BB") {
				if ($new eq "CC") {
					$mark2++;
				}
			}
		}
	}
	@seq=@nextseq;
	@nextseq=qw();
	#############################################################

	########################selection##############################
	my %sel;
	my $sumf=0;
	my $n2=1;
	foreach my $key (@seq) {
		$n2++;
		$key=~m/(.)(.)/;
		my $p1=$1;
		my $p2=$2;
		my $f=1;
		if ($p1 eq "C" and $p2 eq "C") {
		}
		else {
			if ($p1 ne "A") {
				$f=$f*$fitp;
			}
			if ($p2 ne "A") {
				$f=$f*$fitp;
			}
		}
		my $tmp=$key.$n2;
		$sel{$tmp}=$f;
		$sumf=$sumf+$f;
	}
	my @living=qw();
	foreach my $key (keys %sel) {
		my $po=$sel{$key}*$individual/$sumf;
		my $j=rand(1);
		$key=~s/\d//g;
		if ($j<=$po) {
			push @living,$key;
		}
	}
	@seq=@living;
	##############################################################

	###########################record events###############################
	my $nc=0;
	my $na=0;
	my $all=@seq;
	foreach my $k (@seq) {
		if ($k eq "AA") {
			$na++
		}
		if ($k eq "CC") {
			$nc++
		}
	}
	my $ra=$na/$all*100;
	my $rc=$nc/$all*100;
	$probA[$n-1]=$ra;
	$probC[$n-1]=$rc;

	$n++;
	##############################################################
}
my $minp=min(@probA);
my $maxp=max(@probC);

print "$mark1\t$mark2\t$mark";