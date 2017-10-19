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


if ($#ARGV<3) {
    print "USE: smooth_cube.pl INPUT1.CUBE.FITS NX NY OUTPUT.CUBE.FITS\n";
    exit;
}

$infile=$ARGV[0];
$NX=$ARGV[1];
$NY=$ARGV[2];
$iNY=int($NY);
$iNX=int($NX);
$outfile=$ARGV[3];
$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny,$nz)=$a_in->dims;
$kernel=zeroes($iNX,$iNY);
$sum=0;
for ($i=0;$i<$iNX;$i++) {
    for ($j=0;$j<$iNY;$j++) {
	$r=sqrt(($i-$NX/2)**2+($j-$NY/2)**2);
	$val=exp(-0.5*(($r/($NX+$NY)/(2.5*2.345)))**2);
	$sum=$sum+$val;
	set($kernel,$i,$j,$val);
    }
}
$a_out=zeroes($nx,$ny,$nz);
#$kernel=$kernel/$sum;

for ($k=0;$k<$nz;$k++) {
    my $a_in_map=$a_in->slice(":,:,$k");
    my $a_out_map=$a_out->slice(":,:,$k");
    $a_out_map .= med2d($a_in_map,$kernel,{Boundary => Default});    
}


$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
