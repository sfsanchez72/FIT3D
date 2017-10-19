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

if ($#ARGV<4) {
    print "USE: dust_det.pl Rv ratio_obs_21 ratio_expected_21 lambda_1 lambda_2\n";
    exit;
}

$Rv=$ARGV[0];
$rat_obs=$ARGV[1];
$rat_exp=$ARGV[2];
$lambda1=$ARGV[3];
$lambda2=$ARGV[4];

$a_1=A_l($Rv,$lambda1);
$a_2=A_l($Rv,$lambda2);

$Av=(2.5*log10($rat_obs/$rat_exp))/($a_1-$a_2);

#print "$Av $a_1 $a_2\n";
print "$Av\n";

exit;
