#!/usr/bin/perl
#
#
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

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<4) {
    print "USE: correct_map.pl Rv Flux.fits lambda Avmap.fits OUT.fits\n";
    exit;
}

$Rv=$ARGV[0];
$Fluxfile=$ARGV[1];
$lambda=$ARGV[2];
$Avfile=$ARGV[3];
$OUTfile=$ARGV[4];
$Flux=rfits($Fluxfile);
$Av=rfits($Avfile);
$Al=A_l($Rv,$lambda);
$Al=$Al*$Av;
$fac=10**(0.4*$Al);
($nx,$ny)=$Av->dims;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$Av->at($i,$j);
	if ($val<0.001) {
	    set($fac,$i,$j,1);
	}
	$val2=$fac->at($i,$j);
    }
}
$OUT=$Flux*$fac;
$OUT->wfits($OUTfile);
#print "$Flux2\n";
exit;
