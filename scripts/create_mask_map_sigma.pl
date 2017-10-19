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

if ($#ARGV<3) {
    print "USE: create_mask_map_sigma.pl flux_map.fits sigma_map.fits Nsigma mask.fits\n";
    exit;
}

$input_file=$ARGV[0];
$input_sfile=$ARGV[1];
$cut=$ARGV[2];
$output_file=$ARGV[3];

$pdl_input=rfits($input_file);
$pdl_sinput=rfits($input_sfile);
$h=$pdl_input->gethdr;
($nx,$ny)=$pdl_input->dims;
$pdl_output=zeroes($nx,$ny);
$pdl_output->sethdr($h);

for ($j=0;$j<$ny;$j++){
    for ($i=0;$i<$nx;$i++){
	$val=$pdl_input->at($i,$j);
	$sigma=$pdl_sinput->at($i,$j);
	if ($val!=0) {
	    $rat=$val/$sigma;
	    if ($rat>$cut) {
		set($pdl_output,$i,$j,1);
	    }
	} else {
	    if ($sigma<0.0001) {
		set($pdl_output,$i,$j,1);
	    }
	}
    }
}

$pdl_output->wfits($output_file);
exit;


