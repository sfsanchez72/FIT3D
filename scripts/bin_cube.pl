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
    print "USE: bin_cube.pl INPUT.CUBE.fits BIN_BOX OUTPUT.CUBE.fits \n";
    exit;
}

$infile=$ARGV[0];
$N=$ARGV[1];
$outfile=$ARGV[2];
$pdl=rfits($infile);
$h=$pdl->hdr;
($nx,$ny,$nz)=$pdl->dims;
$NX=int($nx/$N);
$NY=int($ny/$N);
$out=zeroes($NX,$NY,$nz);

$$h{NAXIS1}=$NX;
$$h{NAXIS2}=$NY;
$$h{CRPIX1}=$$h{CRPIX1}/$N;
$$h{CRPIX2}=$$h{CRPIX2}/$N;
$$h{CDELT1}=$$h{CDELT1}*$N;
$$h{CDELT2}=$$h{CDELT2}*$N;
$$h{CD1_1}=$$h{CD1_1}*$N;
$$h{CD2_2}=$$h{CD2_2}*$N;



for ($i=0;$i<$NX;$i++) {
    for ($j=0;$j<$NY;$j++) {
	$nx1=$i*$N;
	$nx2=($i+1)*$N;
	$ny1=$j*$N;
	$ny2=($j+1)*$N;
	if ($nx1<0) {
	    $nx1=0;
	}
	if ($ny1<0) {
	    $ny1=0;
	}
	if ($nx2>=$nx) {
	    $nx2=$nx-1;
	}
	if ($ny2>=$ny) {
	    $ny2=$ny-1;
	}

	my $sec=$pdl->slice("$nx1:$nx2,$ny1:$ny2,");
	my $sum=sumover($sec);
	my $sum2=sumover($sum);
	$t = $out->slice("($i),($j),");
	$t .=$sum2/($N*$N);
    }
}

$out->sethdr($h);
$out->wfits($outfile);

exit;
