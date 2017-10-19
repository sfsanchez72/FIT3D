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
if ($#ARGV<1) {
    print "USE: add_header.pl FILE.FITS HEADER_FILE\n";
    print "Header file is a pipe separated file including the keyword and its value\n";
    exit;
}

$file=$ARGV[0];
$hdrfile=$ARGV[1];


$pdl=rfits($file);
open(FH,"<$hdrfile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	($header,$value)=split(/\|/,$line);
	$pdl->hdr->{$header}=$value;
    }
}
close(FH);
$pdl->wfits($file);

exit;

