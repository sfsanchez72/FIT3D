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
use PDL::Transform;


#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<2) {
    print "USE: scale_img.pl INPUT.IMG.fits SCALE_FACTOR OUTPUT.IMG.fits \n";
    exit;
}

$infile=$ARGV[0];
$N=$ARGV[1];
$outfile=$ARGV[2];
$pdl=rfits($infile);
$h=$pdl->hdr;
($nx,$ny)=$pdl->dims;
$NX=int($nx*$N);
$NY=int($ny*$N);
$$h{NAXIS1}=$NX;
$$h{NAXIS2}=$NY;
$$h{CDELT1}=$N*$$h{CDELT1};
$$h{CDELT2}=$N*$$h{CDELT2};
$$h{CRVAL1}=$N*$$h{CRVAL1};
$$h{CRVAL2}=$N*$$h{CRVAL2};

$out=zeroes($NX,$NY);
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$pdl->at($i,$j);
	for ($ii=$i*$N;$ii<($i+1)*$N;$ii++) {
	    for ($jj=$j*$N;$jj<($j+1)*$N;$jj++) {
		set($out,$ii,$jj,$val);
	    }

	}
    }
}

#system("rm -f $outfile");
$out->sethdr($h);
$out->wfits($outfile);

exit;
