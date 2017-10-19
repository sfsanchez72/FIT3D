#!/usr/bin/perl

use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;

#use POSIX;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<4) {
    print "USE:get_index_montecarlo.pl SPEC.txt noise.txt NSIM REDSHIFT DEV\n";    
    exit;
}
$infile=$ARGV[0];
$noisefile=$ARGV[1];
$NSIM=$ARGV[2];
$z=$ARGV[3];
$dev=$ARGV[4];


$spec_file="montecarlo.rss.fits";


open(FH,"<$infile");
open(FH2,"<$noisefile");
while($line=<FH>) {
    chop($line);
    $line2=<FH2>;
    chop($line2);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$id[$n]=$data[0];
	$w[$n]=$data[1];
	$f[$n]=$data[2];
	@data2=split(" ",$line2);
	$noise[$n]=$data[2];
	if ($min>$f[$n]) {
	    $min=$f[$n];
	}
	if ($max<$f[$n]) {
	    $max=$f[$n];
	}
	$n++;
    }
}
close(FH2);
close(FH);


$cdelt=$w[1]-$w[0];
$fwhm=2.345*3*$cdelt;
if ($fwhm<5) {
    $fwhm=5;
}

$call="montecarlo.pl ".$infile." ".$noisefile." ".$NSIM." ".$fwhm." montecarlo.rss.fits";
system($call);

$call="get_index_rss.pl montecarlo.rss.fits ".$z." ".$dev." > get_index_montecarlo.out";
system($call);

$nn=0;
open(FH,"<get_index_montecarlo.out");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$ni=($#data)/2;
	for ($i=0;$i<$ni;$i++) {
	    $val[$i][$nn]=$data[2*$i+2];
	    $name[$i]=$data[2*$i+1];
	}
	$nn++;
    }
}
close(FH);

for ($i=0;$i<$ni;$i++) {
    my @data;
    for ($j=0;$j<$nn;$j++) {
	$data[$j]=$val[$i][$j];
    }
    $med[$i]=median(@data);
    $sig[$i]=sigma(@data);
    print "$name[$i] $med[$i] $sig[$i] ";
}

#print "NN=$n\n";
$med_w=median(@w);
$length=($w[$n-1]-$w[0])/20;
$k=0;
for ($i=0;$i<$n;$i++) {
    if (($w[$i]>($med_w-$length))&&($w[$i]<($med_w+$length))) {
	$S[$k]=$f[$i];
	$N[$k]=$noise[$i];
	$k++;
    }
}
$med_S=median(@S);
$sig_N=sigma(@N);

print "SN  $med_S  $sig_N\n";

exit;
