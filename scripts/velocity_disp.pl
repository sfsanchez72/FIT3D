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


if ($#ARGV<4) {
    print "USE: velocity_disp.pl INPUT_disp_AA.FITS INST_PROF_FWHM WAVELENGTH VELOCITY.fits OUTPUT.FITS\n";
    exit;
}

$infile=$ARGV[0];
$inst=$ARGV[1];
$wave=$ARGV[2];
$vel=$ARGV[3];
$outfile=$ARGV[4];

$a_in=rfits($infile);
$vel_in=rfits($vel);
$wave_in=$wave*(1+$vel_in/300000);
$h=$a_in->gethdr;
$a_out=(sqrt(abs($a_in**2-$inst**2))/$wave_in)*300000;

($nx,$ny)=$a_in->dims;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$a_in->at($i,$j);
	if ($val<$inst) {
	    set($a_out,$i,$j,0);
	}
    }
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
