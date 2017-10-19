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


if ($#ARGV<0) {
    print "USE: img2spectra.pl IMAGE.FITS\n";
    exit;
}
$image=$ARGV[0];

$pdl=rfits($image);
($nx,$ny)=$pdl->dims;

for ($i=0;$i<$ny;$i++) {
    $name="NAME".$i;
    $file=$pdl->hdr->{$name};
    $call="img2spec.pl ".$image." ".$i." ".$file;
    system($call);
}
exit;







