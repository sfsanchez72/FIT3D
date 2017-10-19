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


if ($#ARGV<3) {
    print "USE: smooth_ff.pl INPUT1.FITS NX NY OUTPUT.FITS\n";
    exit;
}

$infile=$ARGV[0];
$NX=$ARGV[1];
$NY=$ARGV[2];
$outfile=$ARGV[3];
$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
$kernel=zeroes($NX,$NY);
$sum=0;
for ($i=0;$i<$NX;$i++) {
    for ($j=0;$j<$NY;$j++) {
	$r=sqrt(($i-$NX/2)**2+($j-$NY/2)**2);
	$val=exp(-0.5*(($r/($NX+$NY)/2))**2);
	$sum=$sum+$val;
	set($kernel,$i,$j,$val);
    }
}
#$kernel=$kernel/$sum;

#$a_out=med2df($a_in,$NX,$NY,{Boundary => Reflect});
$a_out=med2d($a_in,$kernel,{Boundary => Default});
$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
