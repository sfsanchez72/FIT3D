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

if ($#ARGV<4) {
    print "USE: fit_2D_imag.pl INPUT_FILE.FITS MASK.fits[0=MASK] CONFIG_FILE OUTPUT_MOD.FITS DEVICE [min max]\n";
    exit;
}

$input=$ARGV[0];
$mask=$ARGV[1];
$config=$ARGV[2];
$output=$ARGV[3];
$dev=$ARGV[4];
if ($#ARGV==6) {
    $min=$ARGV[5];
    $max=$ARGV[6];
}

$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
$pdl_mask=rfits("$mask");
#($nx,$ny)=$pdl->dims;

$k=1;
open(FH,">slice.junk");
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$pdl->at($i,$j);
	if ($val eq "nan") {
	    $val=0;
	}
	$mask_val=$pdl_mask->at($i,$j);
	if ($mask_val!=0) {
	    print FH "$k $i $j $val\n";
	    $k++;
	}
    }
}
close(FH);

$call="fit_2D_map slice.junk ".$config;
system($call);

$call="plot_out_fit_2D_map.pl slice.junk fit_2D_map.mod ".$dev." R 1";
if ($#ARGV==5) {
    $call=$call." -1 ".$min." ".$max;
}
system($call);

$call="slicegrid_grid.pl fit_2D_map.mod ".$output." ".$nx." ".$ny." 1";
system($call);

$call="slicegrid_grid.pl fit_2D_map.res res_".$output." ".$nx." ".$ny." 1";
system($call);


exit;
