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

if ($#ARGV<1) {
    print "USE: img2slice.pl INPUT_FILE.FITS OUTPUT_SLICE.txt\n";
    exit;
}

$input=$ARGV[0];
$output=$ARGV[1];

$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;

$k=1;
open(FH,"$output");
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$pdl->at($i,$j);
	print FH "$k $i $j\n";
	$k++;
    }
}
close(FH);



exit;
