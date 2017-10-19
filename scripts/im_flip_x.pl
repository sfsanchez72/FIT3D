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
    print "USE: im_flip_x.pl INPUT1.FITS OUTPUT.FITS\n";
    exit;
}

$infile1=$ARGV[0];
$outfile=$ARGV[1];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;
($nx,$ny)=$a_in1->dims;
$nx=$nx-1;
$a_in1=$a_in1->slice("$nx:0,");
$a_in1->sethdr($h);
$a_in1->wfits($outfile);

exit;
