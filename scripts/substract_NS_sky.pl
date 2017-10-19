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
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<2) {
    print "USE: substract_NS_sky.pl INPUT.FITS OUTPUT.FITS CHARGE_OFFSET\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$charge_offset=$ARGV[2];

$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
if ($ny!=(2*int($ny/2))) {
    print "$ny must be even\n";
    exit;
}
$a_out=zeroes($nx,$ny);
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val_obj=$a_in->at($i,$j);
	$j_sky=$j+$charge_offset;
	if ($j_sky<0) {
	    $j_sky=0;
	}
	if ($j_sky>=$ny) {
	    $j_sky=$ny-1;
	}
	$val_sky=$a_in->at($i,$j_sky);
	$val=$val_obj-$val_sky;
#	$k=$j/2;
	set($a_out,$i,$j,$val);
    }
    print "$j/$ny\n";
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
