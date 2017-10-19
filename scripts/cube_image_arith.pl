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


if ($#ARGV<2) {
    print "USE: recal_cube.pl CUBE.fits OPERATOR IMAGE.FITS OUTPUTCUBE.fits\n";
    exit;
}

$infile=$ARGV[0];
$operator=$ARGV[1];
$image=$ARGV[2];
$outfile=$ARGV[3];


$pdl=rfits($infile);
($nx,$ny,$nz)=$pdl->dims;
$h=$pdl->gethdr;

$pdl_image=rfits($image);
($NX,$NY)=$pdl_image->dims;

if (($NX!=$nx)&&($NY!=$ny)) {
    print "Dimensions does not match ($NX,$NY)!=($nx,$ny)\n";
    exit;
}

for ($i=0;$i<$nz;$i++) {
    my $slice=$pdl->slice(",,($i)");
    if ($operator eq "+") {
	$slice.=$slice+$image;
    }

    if ($operator eq "-") {
	$slice.=$slice-$image;
	
    }

    if ($operator eq "/") {
	$slice.=$slice/$image;
	
    }
    if ($operator eq "*") {
	$slice.=$slice*$image;
    }
}



$pdl->wfits($outfile);

exit;
