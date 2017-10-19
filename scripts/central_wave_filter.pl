#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
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

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{'PGPLOT_ENVOPT'}="V";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
if ($#ARGV<0) {
    print "USE: mag_filter.pl filter.dat \n";
    exit;
}

$C=299792.458;

$filter=$ARGV[0];


$wc=0;
$sum=0;
open(FH,"<$filter");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
	$flux_f=$data[2];
	$wave_f=$data[1];
    } else {
	$flux_f=$data[1];
	$wave_f=$data[0];
    }
    $wc=$wc+$wave_f*$flux_f;
    $sum=$sum+$flux_f;
}
close(FH);

$wc=$wc/$sum;

print "$wc\n";

exit;
