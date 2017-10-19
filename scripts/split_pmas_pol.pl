#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module

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



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<2) {
    print "USE: split_pmas_pol.pl INPUT POL1.fits POL2.fits\n";
    exit;
}

$inputfile=$ARGV[0];
$pol1file=$ARGV[1];
$pol2file=$ARGV[2];

$pdl=rfits($inputfile);
($nx,$ny)=$pdl->dims;
$h=$pdl->gethdr;
$NY=$ny/2;
$pdl1=zeroes($nx,$NY);
$pdl2=zeroes($nx,$NY);

for ($j=0;$j<$NY;$j++) {
    $j1=2*$j;
    $j2=2*$j+1;
    for ($i=0;$i<$nx;$i++) {
	$val1=$pdl->at($i,$j1);
	$val2=$pdl->at($i,$j2);
	set($pdl1,$i,$j,$val1);
	set($pdl2,$i,$j,$val2);
    }
}
$$h{NAXIS2}=$NY;
$pdl1->sethdr($h);
$pdl2->sethdr($h);
$pdl1->wfits($pol1file);
$pdl2->wfits($pol2file);


exit;
