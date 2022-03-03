#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use 5.010;
use File::Basename;
use File::Path;
use List::Util qw(min max);
use Getopt::Long;
use Statistics::Descriptive;
use Statistics::Distributions;

#perl do.pl -r 0.1 -f 0.5 -a 1.1
my ($rate,$fitp,$fita);#mutation rate, preserved fitness, fitness advantage
GetOptions(
        'r|rate=f'  => \$rate,
        'f|fitp=f'  => \$fitp,
		'a|fita=f'  => \$fita,
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
my @sample=qw();
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
			$f=$f*$fita;
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


	###########################sampling###############################
	@living=List::Util::shuffle @living;
	@sample=(@sample,@living[0..4]);
	##############################################################
	$n++;
}


@seq=List::Util::shuffle @sample;
$n=0;
my $c1="";
my $c2="";
my @ajs=qw(D E F G H I J K L M N O P Q R S T U);
while ($n<100) {
	$seq[$n]=~m/(.)(.)/;
	my $t1=$1;
	my $t2=$2;
	if ($t1 eq "B") {
		$t1=$ajs[int(rand(18))];
	}
	if ($t2 eq "B") {
		$t2=$ajs[int(rand(18))];
	}
	$c1=$c1.$t1;
	$c2=$c2.$t2;
	$n++;
}
my $r=&get_gft($c1,$c2);#change according to purpose
my $hxy=get_hxy($c1,$c2);#change according to purpose

print "$hxy\t$r";

#################get_hxy############
sub get_hxy{
	my $seq1;
	my $seq2;
	($seq1,$seq2)=@_;
	my @ajs1=split //,$seq1;
	my @ajs2=split //,$seq2;
	my %co1;
	my %co2;
	my $hx=0;
	my $hy=0;
	my $logn=log(20);
	my $len=length ($seq1);
	foreach my $a1 (@ajs1) {
		$co1{$a1}++;
	}
	foreach my $a2 (@ajs2) {
		$co2{$a2}++;
	}
	foreach my $a1 (keys %co1) {
		$hx=$hx-$co1{$a1}/$len*log($co1{$a1}/$len)/$logn;
	}
	foreach my $a2 (keys %co2) {
		$hy=$hy-$co2{$a2}/$len*log($co2{$a2}/$len)/$logn;
	}
	my $hxy=($hx*$hy)**0.5;
	$hxy
}
#################################################


#################get_mi############
sub get_mi{
	my $seq1;
	my $seq2;
	($seq1,$seq2)=@_;
	my @ajs1=split //,$seq1;
	my @ajs2=split //,$seq2;
	my %co1;
	my %co2;
	my $hx=0;
	my $hy=0;
	my $hxy=0;
	my $logn=log(20);
	my $len=length ($seq1);
	foreach my $a1 (@ajs1) {
		$co1{$a1}++;
	}
	foreach my $a2 (@ajs2) {
		$co2{$a2}++;
	}
	foreach my $a1 (keys %co1) {
		$hx=$hx-$co1{$a1}/$len*log($co1{$a1}/$len)/$logn;
	}
	foreach my $a2 (keys %co2) {
		$hy=$hy-$co2{$a2}/$len*log($co2{$a2}/$len)/$logn;
	}
	
	my $n=0;
	my %pair;
	while ($ajs1[$n] and $ajs2[$n]) {
		$pair{"$ajs1[$n],$ajs2[$n]"}++;
		$n++
	}
	foreach my $c (keys %pair) {
		my $pxy=$pair{$c}/$len;
		$hxy=$hxy-$pxy*log($pxy)/$logn;
	}
	my $mi=$hx+$hy-$hxy;
	$mi;
}
#################################################

#################get_mc############
sub get_mc{
	my $seqA;
	my $seqB;
	($seqA,$seqB)=@_;

	#######################################################################
	my @siteA=split //,$seqA;
	my @siteB=split //,$seqB;

	my @disA;
	foreach my $a (@siteA) {
		foreach my $b (@siteA) {
			my $t;
			if ($a eq $b) {
				$t=1;
			}
			else {
				$t=$fitp;
			}
			push (@disA,($t));
		}
	}

	my @disB;
	foreach my $a (@siteB) {
		foreach my $b (@siteB) {
			my $t;
			if ($a eq $b) {
				$t=1;
			}
			else {
				$t=$fitp;
			}
			push (@disB,($t));
		}
	}

	my $dataA = Statistics::Descriptive::Full->new();
	my $dataB = Statistics::Descriptive::Full->new();
	$dataA->add_data(\@disA);
	$dataB->add_data(\@disB);

	my $meanA=$dataA->mean();
	my $meanB=$dataB->mean();
	my $sdA=$dataA->standard_deviation();
	my $sdB=$dataB->standard_deviation();

	my $n=0;
	my $sum=0;
	my $r=0;
	if ($sdA ne 0 and $sdB ne 0) {
		foreach my $v (@disA) {
			my $tmp=($disA[$n]-$meanA)*($disB[$n]-$meanB)/($sdA*$sdB);
			$sum=$sum+$tmp;
			$n++;
		}
		$r=$sum/$n;
	}
	#######################################################################

	$r;
}
#################################################


#################get_gft############
sub get_gft{
	my $seq1;
	my $seq2;
	($seq1,$seq2)=@_;
	my @ajs1=split //,$seq1;
	my @ajs2=split //,$seq2;

	my %co1;
	my %co2;
	my $len=length ($seq1);
	foreach my $a1 (@ajs1) {
		$co1{$a1}++;
	}
	foreach my $a2 (@ajs2) {
		$co2{$a2}++;
	}
	my %co12;
	my $n=0;
	foreach my $a (@ajs1) {
		my $b=$ajs2[$n];
		my $ab=$a.$b;
		$co12{$ab}++;
		$n++;
	}

	my $sum=0;
	my $df=0;
	foreach my $v1 (keys %co1) {
		foreach my $v2 (keys %co2) {
			my $T=($co1{$v1}*$co2{$v2})/$len;
			my $O;
			my $v=$v1.$v2;
			if ($co12{$v}) {
				$O=$co12{$v};
			}
			else {
				$O=0;
			}
			my $tmp=($O-$T)**2/$T;
			$sum=$sum+$tmp;
			$df++;
		}
	}


	my $x2=$sum;
	$df=$df-1;
	my $p=0;
	if ($df ne 0) {
		$p=1-Statistics::Distributions::chisqrprob ($df,$x2);
	}
	$p;
}
#################################################