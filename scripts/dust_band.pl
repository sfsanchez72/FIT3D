#!/usr/bin/perl
#
#
#
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

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<2) {
    print "USE: dust_band.pl A_band_1 lambda_1  lambda_2\n";
    exit;
}

$A1=$ARGV[0];
$lambda1=$ARGV[1];
$lambda2=$ARGV[2];

$a_1=A_l2(4.05,$lambda1);
$a_2=A_l2(4.05,$lambda2);

#print "$a_2 $a_1 ";
$A2=$A1*$a_2/$a_1;

print "$A2\n";

exit;
