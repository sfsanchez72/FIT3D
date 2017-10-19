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


if ($#ARGV<1) {
    print "USE: imarith.pl INPUT.FITS LOG10_INPUT.fits\n";
    exit;
}

$infile1=$ARGV[0];
$outfile=$ARGV[1];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;

($nx,$ny)=$a_in1->dims;
$pdl=zeroes($nx,$ny);
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$a_in1->at($i,$j);
	if ($val>0) {
	    $l=log10($val);
	    set($pdl,$i,$j,$l);
	} else {
	    set($pdl,$i,$j,-1000);
	}
    }
}

$pdl->sethdr($h);
$pdl->wfits($outfile);

exit;
