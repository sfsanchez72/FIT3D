#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
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
use PDL::Slatec;
use PDL::Image2D;

use PDL::NiceSlice;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<2) {
    print "USE: create_sky_vimos_med.pl INPUT mask_file OUTPUT [NBIN]\n";
    exit;
}

$inputfile=$ARGV[0];
$maskfile=$ARGV[1];
$outputfile=$ARGV[2];
$nbin=400;
if ($#ARGV==3) {
    $nbin=$ARGV[3];
}

$n=0;
open(FH,"<$maskfile");
while($line=<FH>) {
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $mask[$n]=$data[1];
    $n++;
}
close(FH);

$input=rfits($inputfile);
($n1,$n2)=$input->dims;

if ($n2!=$n) {
    print "The Y-axis of $inputfile is different than the number\n of lines in the mask file\n";
    exit;
}


$nk=$n/$nbin;
#$r_input=$input->xchg(0,1);
#$r_output=$input->xchg(0,1);
#$sky=zeroes($n1,$nk);
for ($k=0;$k<$nk;$k++) {
    $nsky=0;
#	my @sky;
    for ($j=($k*$nbin);$j<(($k+1)*$nbin);$j++) {
	if ($mask[$j]==1) {
	    $nsky++;
	}
    }
    print "$nsky/$nbin\n";
    if ($nsky>0) {
	$sky=zeroes($n1,$nsky);
	$nsky=0;
	for ($j=($k*$nbin);$j<(($k+1)*$nbin);$j++) {
	    if ($mask[$j]==1) {
		$t=$sky->slice(":,$nsky"); $t .=$t+($input->slice(":,$j"));
		$nsky++;
	    } 
	}
#	$t=$sky->slice(":,$k"); $t .=$t/$nsky;

	$spectrum=medover($sky->xchg(0,1));
    #$spectrum=medover($sky);
	@dims=$sky->dims();
	@dims_spec=$spectrum->dims();
#    print "@dims @dims_spec\n";
    }
    for ($j=($k*$nbin);$j<(($k+1)*$nbin);$j++) {
	if ($mask[$j]==1) {
	    $t=$input->slice(":,$j"); $t .=$t-$spectrum;
	} else {
	    $t=$input->slice(":,$j"); $t .=0;
	}
    }    
}
$input->wfits($outputfile);

exit;
