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
    print "USE: cr_mean_clean.pl input.fits mean.fits X_BOX Y_BOX lower_limit OUTPUT.FITS\n";
    exit;
}

$infile=$ARGV[0];
$meanfile=$ARGV[1];
$x_box=$ARGV[2];
$y_box=$ARGV[3];
$limit=$ARGV[4];
$outfile=$ARGV[5];

$pdl=rfits($infile);
$h=$pdl->gethdr;
($nx,$ny)=$pdl->dims();
$pdl_mean=rfits($meanfile);
$pdl_sub=$pdl-$pdl_mean;

$smoothed=med2d($pdl_sub, ones($x_box,$y_box), {Boundary => Reflect});

$pdl_comp=$pdl_sub-$smoothed;

$pdl_out=$pdl;
for ($j=0;$j<$ny;$j++) {
    for($i=0;$i<$nx;$i++) {
	$val=$pdl_comp->at($i,$j);
	if ($val>$limit) {	    
	    $mean=$smoothed->at($i,$j);
	    set($pdl_out,$i,$j,$mean);
	}
    }
}

$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);

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
