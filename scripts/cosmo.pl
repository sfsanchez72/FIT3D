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

if ($#ARGV<0) {
    print "USE: cosmo.pl redshift\n";
    exit;
}


$redshift=$ARGV[0];

@data=cosm($redshift,71,0.27,0.73);

print "$data[5] Kpc/arcsec\n";

$m_abs=m_abs(71,0.5,$redshift,0,0);
#print "m_abs=m+$m_abs\n";
#$m_abs=m_abs_low(71,0.5,$redshift,0,0);
#print "m_abs=m+$m_abs\n";
#$m_abs=m_abs_lambda(71,0.3,0.7,$redshift,0,0);
#print "m_abs=m+$m_abs\n";
exit;
