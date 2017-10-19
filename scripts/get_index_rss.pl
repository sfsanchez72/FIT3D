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

#use POSIX;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<2) {
    print "USE:get_index_rss.pl SPEC.RSS.fits REDSHIFT DEV\n";    
    exit;
}

$spec_file=$ARGV[0];
$z=$ARGV[1];
$dev=$ARGV[2];

$a=rfits($spec_file);
($nx,$ny)=$a->dims;
for ($j=0;$j<$ny;$j++) {
    $call="img2spec.pl ".$spec_file." ".$j." get_index_rss.spec";
    system($call);
    $call="get_index.pl get_index_rss.spec ".$z." ".$dev;
    system($call);
}


exit;
