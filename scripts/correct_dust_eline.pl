#!/usr/bin/perl
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
    print "USE: correct_line.pl Flux lambda Av\n";
    exit;
}

$Flux=$ARGV[0];
$lambda=$ARGV[1];
$Av=$ARGV[2];
#$Al=A_l(5.5,$lambda);
$Al=A_l(3.1,$lambda);
#print "$Al ";
$Al=$Al*$Av;
$fac=10**(0.4*$Al);
$Flux2=$Flux*$fac;
#print "$Al $fac $Flux2\n";
print "$Flux2\n";
exit;
