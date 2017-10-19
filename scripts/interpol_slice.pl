#!/usr/bin/perl
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
    print "USE: interpol_slice.pl slice.txt image_out.fits dpix grid_func grid_option\n";
    exit;
}

$slice=$ARGV[0];
$out_img=$ARGV[1];
$dpix=$ARGV[2];
$gf=$ARGV[3];
$go=$ARGV[4];

open(FH,"<$slice");
open(OUT,">slice.junk");
$line=<FH>;
while($line=<FH>) {
    print OUT "$line";
}
close(OUT);
close(FH);

$call="interpol -if slice.junk -of ".$out_img." -dp ".$dpix." -gf ".$gf." -go ".$go;
system($call);

exit;



