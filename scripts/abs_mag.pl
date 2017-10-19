#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
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

if ($#ARGV<1) {
    print "USE: abs_mag.pl redshift magnitude [Evolution Passive]\n";
    exit;
}


$redshift=$ARGV[0];
$mag=$ARGV[1];

$out=m_abs_lambda(71,0.27,0.73,$redshift,$mag,0);

if ($#ARGV==2) {
    $out=$out+4*log10(1+$redshift);
}

print "m_abs = $out mag\n";


exit;
