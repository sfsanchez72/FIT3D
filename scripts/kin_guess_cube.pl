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

if ($#ARGV<4) {
    print "USE: kin_guess_cube.pl cube.fits start_wavelength end_wavelength out_vel_file DEVICE\n";
    exit;
}

$input_cube=$ARGV[0];
$start_w=$ARGV[1];
$end_w=$ARGV[2];
$out_file=$ARGV[3];
$device=$ARGV[4];

$pdl_input_cube=rfits($input_cube);
$h=$pdl_input_cube->gethdr;
$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$ny=$pdl_input_cube->hdr->{"NAXIS2"};
$nz=$pdl_input_cube->hdr->{"NAXIS3"};
$crpix=$pdl_input_cube->hdr->{"CRPIX3"};
$crval=$pdl_input_cube->hdr->{"CRVAL3"};
$cdelt=$pdl_input_cube->hdr->{"CDELT3"};
$pdl_input_cube->inplace->setnantobad; 
$pdl_input_cube->inplace->setbadtoval(-1); 

$map_vel=zeroes($nx,$ny);
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	my @W;
	my @F;
	my @Fcen;
	$y_min=1e12;
	$y_max=-1e12;
	for ($k=0;$k<$nz;$k++) {
	    $wave=$crval+$cdelt*($k+1-$crpix);
	    $val=$pdl_input_cube->at($i,$j,$k);
	    $W[$k]=$wave;
	    $F[$k]=$val;
	    $Fcen[$k]=$pdl_input_cube->at(int($nx/2),int($ny/2),$k);	    
	    if ($y_min>$F[$k]) {
		$y_min=$F[$k];
	    }
	    if ($y_max<$F[$k]) {
		$y_max=$F[$k];
	    }
#	    print "$k $W[$k] $F[$k] $Fcen[$k] $y_min $y_max\n";
	}

	$median=median(@F);
	$sigma=sigma(@F);
	$y_peak=-1e12;
	$n_peaks=0;
#	my @W_peaks;
	
	for ($k=1;$k<$nz-1;$k++) {
	    if (($F[$k]>($median+$sigma))&&($y_peak<$F[$k])) {
		
	    }
	}

	
	pgbegin(0,$dev,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($start_w,$end_w,$y_min-0.01,$y_max+0.01,0,0);
	pglabel("Wavelength (\\A)","Flux (units)","($i,$j) ($nx,$ny)");
	pgsci(2);
	pgline($nz,\@W,\@Fcen);
	pgsci(1);
	pgline($nz,\@W,\@F);
	pgclose;
	pgend;
#	print "Press Enter ($i/$nx,$j/$ny)\n"; <stdin>;
    }
}



exit;
