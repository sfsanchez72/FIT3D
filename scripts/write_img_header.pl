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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");
if ($#ARGV<2) {
    print "USE: write_img_header.pl FILE.FITS HEADER VALUE\n";
    exit;
}

$file=$ARGV[0];
$header=$ARGV[1];
$value=$ARGV[2];


$pdl=rfits($file);
$h=$pdl->gethdr;
$h2={$h,$header=>$value};
#$$h{$header}=$value;
#$pdl->sethdr($h2);
$pdl->hdr->{$header}=$value;
$pdl->wfits($file);

exit;

