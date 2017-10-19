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
    print "USE: med2df.pl INPUT.FITS OUTPUT.fits X_WIDTH Y_WIDTH\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$dx=$ARGV[2];
$dy=$ARGV[3];
$pdl=rfits($infile);
$mpdl=med2df($pdl,$dx,$dy,{Boundary=>Reflect});
#$mpdl=med2d($pdl,ones($dx,$dy),{Boundary=>Reflect});
$mpdl->wfits($outfile);

exit;


