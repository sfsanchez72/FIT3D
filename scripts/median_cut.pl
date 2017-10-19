#!/usr/bin/perl
#
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
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<4) {
    print "USE: median_spec.pl RSS.fits NX1 NX2 NSIG_MIN NSIG_MAX \n";
    exit;
}

$infile=$ARGV[0];
$nx1=$ARGV[1];
$nx2=$ARGV[2];
$ns1=$ARGV[3];
$ns2=$ARGV[4];
#$fiber_flat=$ARGV[2];


$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
$new_ny=$ny;
if ($#ARGV==2) {
    $new_ny=$ARGV[2];
}


$n=0;
my @cut;

for ($i=$nx1;$i<$nx2;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$cut[$n]=$in_array[$j][$i];
#	print "$n $cut[$n]\n";
	$n++;
    }
}

$m=median(@cut);
$s=sigma(@cut);
$k=0;
for ($i=0;$i<$n;$i++) {
    if (($cut[$i]>($m-$s*$ns1))&&($cut[$i]<($m+$s*$ns2))) {
	$cut2[$k]=$cut[$i];
	$k++;
    }
}


$spec=median(@cut2);



print "$infile $spec $k\n";

exit;
