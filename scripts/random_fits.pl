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


if ($#ARGV<2) {
    print "USE: random_fits.pl INPUT.FITS OUTPUT.fits RANGE\n";
    exit;
}

$infile1=$ARGV[0];
$outfile=$ARGV[1];
$range=$ARGV[2];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;

($nx,$ny)=$a_in1->dims;
$pdl=zeroes($nx,$ny);
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$a_in1->at($i,$j);
	$val=$val-rand($range);
	print "$val\n";
	set($pdl,$i,$j,$val);
    }
}

$pdl->sethdr($h);
$pdl->wfits($outfile);

exit;
