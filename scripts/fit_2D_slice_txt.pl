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

use PDL::Core;
use PDL::Graphics::LUT;
use Carp;

$ENV{'PGPLOT_ENVOPT'}="IV";


#$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
#$galfit="/home/sanchez/sda1/galfit/galfit";
#$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<3) {
    print "USE: fit_2D_imag.pl INPUT_FILE.TXT CONFIG_FILE OUTPUT_MOD.TXT DEVICE [min max] [bias]\n";
    exit;
}

$input=$ARGV[0];
$config=$ARGV[1];
$output=$ARGV[2];
$dev=$ARGV[3];
if ($#ARGV==5) {
    $min=$ARGV[4];
    $max=$ARGV[5];
}

$bias=0;
if ($#ARGV==6) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $bias=$ARGV[6];
}

$y_min=1e12;
$y_max=-1e12;
$n=0;

$call="fit_2D_map ".$input." ".$config." ".$bias;
system($call);

$call="plot_slices_3.pl ".$input." fit_2D_map.mod fit_2D_map.res ".$dev." C 2.68";
if ($#ARGV>4) {
    $call=$call." -1 ".$min." ".$max;
}
system($call);
print "$call\n";


$call="cp fit_2D_map.mod ".$output;
system($call);


exit;
