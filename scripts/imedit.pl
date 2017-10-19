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

if ($#ARGV<5) {
    print "USE: imedit.pl INPUT.FITS OUTPUT.fits lower_limit lower_val upper_limit upper_val\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$lower=$ARGV[2];
$lval=$ARGV[3];
$upper=$ARGV[4];
$uval=$ARGV[5];
$pdl=rfits($infile);
($nx,$ny)=$pdl->dims();
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$pdl->at($i,$j);
	if ($val<$lower) {
	    set($pdl,$i,$j,$lval);
	}
	if ($val>$upper) {
	    set($pdl,$i,$j,$uval);
	}
    }
}

$pdl->wfits($outfile);

exit;


