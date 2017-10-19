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

use PDL::NiceSlice;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<3) {
    print "USE: binning.pl INPUT X_bin Y_bin OUTPUT\n";
    exit;
}

$inputfile=$ARGV[0];
$Nx=$ARGV[1];
$Ny=$ARGV[2];
$outputfile=$ARGV[3];

$pdl_in=rfits($inputfile);
($nx,$ny)=$pdl_in->dims;
$nx_new=int($nx/$Nx);
$ny_new=int($ny/$Ny);
print "$nx_new,$ny_new\n";
$pdl_out=zeroes($nx_new,$ny_new);



for ($i=0;$i<$nx_new;$i++) {
    for ($j=0;$j<$ny_new;$j++) {
	$val=0;
	for ($ii=$i*$Nx;$ii<($i+1)*$Nx;$ii++) {
	    for ($jj=$j*$Ny;$jj<($j+1)*$Ny;$jj++) {
		$val=$val+$pdl_in->at($ii,$jj);
	    }
	}
	set($pdl_out,$i,$j,$val);
    }
}
$pdl_out->wfits($outputfile);


exit;




