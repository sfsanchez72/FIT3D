#!/usr/bin/perl

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

if ($#ARGV<2) {
    print "USE: create_mask_map.pl map.fits cut mask.fits\n";
    exit;
}

$input_file=$ARGV[0];
$cut=$ARGV[1];
$output_file=$ARGV[2];

$pdl_input=rfits($input_file);
$h=$pdl_input->gethdr;
($nx,$ny)=$pdl_input->dims;
$pdl_output=zeroes($nx,$ny);
$pdl_output->sethdr($h);

for ($j=0;$j<$ny;$j++){
    for ($i=0;$i<$nx;$i++){
	$val=$pdl_input->at($i,$j);
	if ($val>$cut) {
	    set($pdl_output,$i,$j,1);
	}
    }
}

$pdl_output->wfits($output_file);
exit;


