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


if ($#ARGV<1) {
    print "USE: imarith.pl INPUT.fits OUTPUT.FITS\n";
    exit;
}

$infile1=$ARGV[0];
$outfile=$ARGV[1];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;
$a_in1=sqrt($a_in1);
$a_in1->sethdr($h);
$a_in1->wfits($outfile);

exit;
