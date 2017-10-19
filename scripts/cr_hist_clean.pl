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


if ($#ARGV<5) {
    print "USE: cr_hist_clean.pl input.fits Nsigma X_BOX Y_BOX clean.fits cr.fits \n";
    exit;
}

$infile=$ARGV[0];
$nsigma=$ARGV[1];
$x_box=$ARGV[2];
$y_box=$ARGV[3];
$outfile=$ARGV[4];
$mask=$ARGV[5];

$pdl=rfits($infile);
$h=$pdl->gethdr;
($nx,$ny)=$pdl->dims();

$pdl_out=$pdl;
$pdl_mask=ones($nx,$ny);
for ($j=$y_box;$j<$ny-$y_box;$j=$j+$y_box) {
    for($i=$x_box;$i<$nx-$x_box;$i=$i+$x_box) {
	$nx1=$i-$x_box;
	$nx2=$i+$x_box;
	$ny1=$j-$y_box;
	$ny2=$j+$y_box;
	if ($nx1<0) {
	    $nx1=0;
	    $nx2=2*$x_box+1;
	}
	if ($ny1<0) {
	    $ny1=0;
	    $ny2=2*$y_box+1;
	}
	if ($nx2>($nx-1)) {
	    $nx2=$nx-1;
	    $nx1=$nx-1-(2*$x_box+1);
	}
	if ($ny2>($ny-1)) {
	    $ny2=$ny-1;
	    $ny1=$ny-1-(2*y_$box+1);
	}
	my $sec=$pdl->slice("$nx1:$nx2,$ny1:$ny2");
	($mean,$rms,$median,$min,$max) = stats($sec);
	my @values=list($sec);
	@sort= sort { $values[$a] <=> $values[$b] } (0..$#values); 
	$found=0;
	$I=0;
	$limit=0;
	for ($ii=1;$ii<$#values+1;$ii++) {
	    $delta=$short[$ii]-$short[$ii-1];
	    if (($short[$ii]>$median)&&($delta>($nsigma*$rms))) {
		$found=1;
		$I=$ii;
		$limit=$short[$ii];
		$ii=1e12;
	    }

	}
	for ($jj=$ny1;$jj<$ny2;$jj++) {
	    for ($ii=$nx1;$ii<$nx2;$ii++) {
		$val=$pdl->at($ii,$jj);
		if ($val>$limit) {
		    set($pdl_mask,$ii,$jj,0);
		}
	    }
	}

	

    }
	print "$j\n";
}

#	    set($pdl_out,$i,$j,$mean);



$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);
$pdl_mask->sethdr($h);
$pdl_mask->wfits($mask);

exit;

sub patata {
    $I1=$i-$x_box;
	$I2=$i+$x_box;
    if ($I1<0) {
	$I2=$i+2*$x_box;
	$I1=0;
    }
    if ($I2>$nx) {
	$I2=$nx;
	$I1=$i-2*$x_box;
    }
    my $sec=$pdl_sub->slice("$I1:$I2,($j)");
    ($mean,$rms,$median,$min,$max) = stats($sec);

}
