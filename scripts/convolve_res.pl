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
    print "USE: convolve_res.pl INPUT.FITS DELTA_RES.FITS OUTPUT.FITS\n";
    exit;
}

$infile1=$ARGV[0];
$infile2=$ARGV[1];
$outfile=$ARGV[2];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;
$a_in2=rfits($infile2);

($nx,$ny)=$a_in1->dims;
$crval=$a_in1->hdr->{CRVAL1};
$cdelt=$a_in1->hdr->{CDELT1};
$pdl_out=zeroes($nx,$ny);

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$sigma=$a_in2->at($i,$j);
	$sigma=$sigma/$cdelt;
	$sum=0;
	$flux=0;
	$imin=$i-int(4*$sigma);
	$imax=$i+int(4*$sigma);
	if ($imin<0) {
	    $imin=0;
	}
	if ($imax>=$nx) {
	    $imax=$nx-1;
	}
	for ($ii=$imin;$ii<$imax;$ii++) {
	    $w=exp(-0.5*(($ii-$i)/$sigma)**2);
	    $sum=$sum+$w;
	    $fin=$a_in1->at($ii,$j);
	    $flux=$flux+$fin*$w;
	}
	if ($sum!=0) {
	    $flux=$flux/$sum;
	} else {
	    $flux=$a_in1->at($i,$j);
	}
#	print "$i,$j $sum\n";
	set($pdl_out,$i,$j,$flux);
    }
}


$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);


exit;
