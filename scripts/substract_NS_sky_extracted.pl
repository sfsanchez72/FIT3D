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


if ($#ARGV<1) {
    print "USE: substract_NS_sky_extracted.pl INPUT.FITS OUTPUT.FITS\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];

$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
if ($ny!=(2*int($ny/2))) {
    print "$ny must be even\n";
    exit;
}
$a_out=zeroes($nx,int($ny/2));
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j=$j+2) {
	$val_obj=$a_in->at($i,$j);
	$val_sky=$a_in->at($i,$j+1);
	$val=$val_obj-$val_sky;
	$k=$j/2;
	set($a_out,$i,$k,$val);
    }
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
