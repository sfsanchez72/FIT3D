#!/usr/bin/perl
#
#
#
use PDL;
use PDL::Core;
use PDL::Graphics::LUT;


use Astro::FITS::CFITSIO;
use PGPLOT;
use Carp;


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



if ($#ARGV<6) {
    print "fits2ppm_hist.pl R.fits G.fits B.fits output.ppm min max gamma\n";
    exit;
}

$Rfile=$ARGV[0];
$Gfile=$ARGV[1];
$Bfile=$ARGV[2];
#$max=$ARGV[3];
$output=$ARGV[3];
$min=$ARGV[4];
$max=$ARGV[5];
$gamma=$ARGV[6];

$pdl_R=rfits($Rfile);
$pdl_G=rfits($Gfile);
$pdl_B=rfits($Bfile);

($nx,$ny)=$pdl_R->dims;

for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val_R=$pdl_R->at($i,$j);
	$val_G=$pdl_G->at($i,$j);	
	$val_B=$pdl_B->at($i,$j);
	if ($val_R!=1*$val_R) {
	    $val_R=1e20;
	}
	if ($val_G!=1*$val_G) {
	    $val_G=1e200;
	}
	if ($val_B!=1*$val_B) {
	    $val_B=1e200;
	}

	if ($val_R eq "BAD") {
	    $val_R=1e200;
	}
	if ($val_G eq "BAD") {
	    $val_G=1e200;
	}
	if ($val_B eq "BAD") {
	    $val_B=1e200;
	}
	if ($val_R>0) {
	    $val_R=log10($val_R);
	} else {
	    $val_R=-1e12;
	}
	if ($val_V>0) {
	    $val_V=log10($val_V);
	} else {
	    $val_V=-1e12;
	}
	if ($val_B>0) {
	    $val_B=log10($val_B);
	} else {
	    $val_B=-1e12;
	}

	set($pdl_R,$i,$j,$val_R);
	set($pdl_G,$i,$j,$val_G);
	set($pdl_B,$i,$j,$val_B);

    }
}


$pdl_R=abs($pdl_R)**$gamma;
$pdl_G=abs($pdl_G)**$gamma;
$pdl_B=abs($pdl_B)**$gamma;


($nx,$ny)=$pdl_R->dims;

$pdl_R=(abs((255/($max-$min))*($pdl_R-$min)));
$pdl_G=(abs((255/($max-$min))*($pdl_G-$min)));
$pdl_B=(abs((255/($max-$min)*1.05)*($pdl_B-$min)));

open(FH,">$output");
print FH "P3\n";
print FH "# $output\n";
print FH "$ny $nx\n";
print "$nx $ny\n";
print FH "256\n";
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$r=int($pdl_R->at($i,$j));
	$g=int($pdl_G->at($i,$j));
	$b=int($pdl_B->at($i,$j));
	if (abs($r)>1e20) {
	    $r=255;
	}
	if (abs($g)>1e20) {
	    $g=255;
	}
	if (abs($b)>1e20) {
	    $b=255;
	}
	if ($r<0) {
	    $r=0;
	}
	if ($g<0) {
	    $g=0;
	}
	if ($b<0) {
	    $b=0;
	}
	if ($r>255) {
	    $r=255;
	}
	if ($g>255) {
	    $g=255;
	}
	if ($b>255) {
	    $b=255;
	}

	if (($r==0)&&($g==0)&&($b==0)) {
#	    $r=255;
#	    $g=255;
#	    $b=255;
	}
	
	print FH "$r $g $b\n";
    }
}
close(FH);
system("display $output");

exit;
