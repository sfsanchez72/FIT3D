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
    print "USE: create_sky_vimos.pl INPUT mask_file OUTPUT\n";
    exit;
}

$inputfile=$ARGV[0];
$maskfile=$ARGV[1];
$outputfile=$ARGV[2];

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


$nk=$n/400;
#$r_input=$input->xchg(0,1);
#$r_output=$input->xchg(0,1);
$sky=zeroes($n1,$nk);
for ($k=0;$k<$nk;$k++) {
    $nsky=0;
    for ($j=($k*400);$j<(($k+1)*400);$j++) {
	if ($mask[$j]==1) {
	    $t=$sky->slice(":,$k"); $t .=$t+($input->slice(":,$j"));
	    $nsky++;
	}
    }
    print "$nsky/400\n";
    $t=$sky->slice(":,$k"); $t .=$t/$nsky;

    for ($j=($k*400);$j<(($k+1)*400);$j++) {
	if ($mask[$j]==1) {
	    $t=$input->slice(":,$j"); $t .=$t-($sky->slice(":,$k"));	    
	} else {
	    $t=$input->slice(":,$j"); $t .=-1000;
	}
    }    
}
$input->wfits($outputfile);

exit;

	    for ($i=0;$i<$n1;$i++){
		$in_val=$sky->at($i,$k);
		$add_val=$input->at($i,$j);
		$in_val=$in_val+$add_val;
		set($sky,$i,$k,$in_val);
		$nsky++;
	    }

	    for ($i=0;$i<$n1;$i++){
		$in_val=$input->at($i,$j);
		$add_val=$sky->at($i,$k);
		$in_val=$in_val-$add_val;
		set($output,$i,$j,$in_val);
	    }
