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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<0) {
    print "USE: get_airmass.pl FITSFILE\n";
    exit;
}
$infile=$ARGV[0];
#@pos=read_img_headers($infile,["RA","DEC"]);
#@pos=read_img_headers($infile,["HIERARCH CAHA TEL POS SET RA","HIERARCH CAHA TEL POS SET DEC"]);
@pos=read_img_headers($infile,["AIRMASS"]);
print "$infile $pos[0]\n";

exit;
