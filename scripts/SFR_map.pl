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

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<3) {
    print "USE: l_abs.pl obs_lum_line.fits factor redshift OUT.FITS [h Om Ol]\n";
    print "lum_line units 10^-16 erg s-1 cm-2\n";
    print "factor = 2.0 for [OII] or 0.89 for Halpha\n";
    exit;
}

$l_in=$ARGV[0];
$l=rfits($l_in);
$factor=$ARGV[1];
$z=$ARGV[2];
$outfile=$ARGV[3];
$h=70.4;
#$h=70.4;
$Om=0.268;
$Ol=0.732;
if ($#ARGV==6) {
    $h=$ARGV[4];
    $Om=$ARGV[5];
    $Ol=$ARGV[6];
}

$COSMO=L_abs_lambda($h,$Om,$Ol,$z,1);
$L=$COSMO*$l;
$SFR=$L*(1e-16)*$factor*(1e-41);

$SFR->wfits($outfile);

exit;
