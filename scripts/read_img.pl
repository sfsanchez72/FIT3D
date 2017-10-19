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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<2) {
    print "USE: read_img.pl INPUT_FILE.FITS X Y \n";
    exit;
}

$input=$ARGV[0];
$x=$ARGV[1];
$y=$ARGV[2];


$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
if (($x>-1)&&($x<$nx)&&($y>-1)&&($y<$ny)) {
    $val=$pdl->at($x,$y);
    print "$val\n";
} else {
    print "$x,$y out of range ($nx,$ny)\n";
}
exit;




