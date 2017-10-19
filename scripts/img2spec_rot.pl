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

if ($#ARGV<1) {
    print "USE: img2spec.pl INPUT_FILE.FITS NX OUTPUTFILE.txt\n";
    exit;
}

$input=$ARGV[0];
$NX=$ARGV[1];
$output=$ARGV[2];


$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};
if ($cdelt==0) {
    $cdelt=1;
}
open(FH,">$output");
for ($i=0;$i<$ny;$i++) {
    $k=$i+1;
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($NX,$i);
    if ($flux[$i] eq "BAD") {
	$flux[$i]=0;
    }
    print FH "$k $wave[$i] $flux[$i]\n";
}
close(FH);

exit;
