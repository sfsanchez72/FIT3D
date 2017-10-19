#!/usr/bin/perl
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
#use PDL::Graphics::TriD;
#use PDL::Graphics::TriD::Image;
use PDL::Fit::Gaussian;
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<3) {
    print "USE: get_spec_cube.pl CUBE.fits X_C Y_C output_spec.txt\n";
    exit;
}

$input=$ARGV[0];
$x_c=$ARGV[1];
$y_c=$ARGV[2];
$output=$ARGV[3];
$pdl_cube=rfits("$input");
($NX,$NY,$NZ)=$pdl_cube->dims();

if (($x_c>=$NX)||($y_c>=$NY)) {
    print "Position out of boundaries\n";
    exit;
}

$h=$pdl_cube->gethdr;
$crval=$pdl_cube->hdr->{"CRVAL3"};
$cdelt=$pdl_cube->hdr->{"CDELT3"};
$crpix=$pdl_cube->hdr->{"CRPIX3"};
open(FH,">$output");
$k=1;
for ($i=0;$i<$NZ;$i++) {
    $wave=$crval+$cdelt*($i+1-$crpix);
    $val=$pdl_cube->at($x_c,$y_c,$i);
    print FH "$k $wave $val\n";
    $k++;
}
close(FH);

exit;
