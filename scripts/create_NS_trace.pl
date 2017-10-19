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


if ($#ARGV<2) {
    print "USE: create_NS_trace.pl INPUT.FITS delta OUTPUT.FITS\n";
    exit;
}

$infile=$ARGV[0];
$delta=$ARGV[1];
$outfile=$ARGV[2];

$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
$a_out=zeroes($nx,2*$ny);
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$a_in->at($i,$j);
	$k=2*$j;
	set($a_out,$i,$k,$val);
	set($a_out,$i,$k+1,$val+$delta);
    }
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
