#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<5) {
    print "USE: sky_sub_cube.pl input_cube.fits Xmin Ymin Xmax Ymax output_cube\n";
    exit;
}


$input_cube=$ARGV[0];
$x0=$ARGV[1]-1;
$y0=$ARGV[2]-1;
$x1=$ARGV[3]-1;
$y1=$ARGV[4]-1;
$output_cube=$ARGV[5];

$pdl_in=rfits($input_cube);
($nx,$ny,$nz)=$pdl_in->dims;
$H=$pdl_in->gethdr;

$nx_0=int($nx*$scale);
$ny_0=int($ny*$scale);

$pdl_out=$pdl_in;

$pdl_sec=$pdl_in->slice("$x0:$x1,$y0:$y1,:");
$spec1=medover($pdl_sec->xchg(0,1));
$spectrum=medover($spec1);



for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$t=$pdl_out->slice("($i),($j),:");
	$t .= $t-$spectrum;
    }
}

$pdl_out->wfits($output_cube);

exit;

