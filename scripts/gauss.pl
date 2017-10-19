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
    print "USE: gauss.pl INPUT.FITS OUTPUT.fits FWHM\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$fwhm=$ARGV[2];
$sig=$fwhm/2.345;
$dx=int($fwhm*3);
if ($dx==(2*int($dx/2))) {
    $dx=$dx+1;
}
$x_c=int($dx/2);
$kernel=zeroes($dx,$dx);
for ($i=0;$i<$dx;$i++) {
    for ($j=0;$j<$dx;$j++) {
	$g=exp(-(($i-$x_c)**2+($j-$x_c)**2)/(2*$sig**2));
	set($kernel,$i,$j,$g);
    }
}
$sum=sum($kernel);
$kernel=$kernel/$sum;
$pdl=rfits($infile);
$mpdl=med2d($pdl,$kernel,{Boundary=>Reflect});

@dat=stats($pdl);
$f1=$dat[0];
@dat=stats($mpdl);
$f2=$dat[0];

$nx=$mpdl->dims;
$mpdl=$mpdl*($f1/$f2);
$h=$pdl->gethdr();
$mpdl->sethdr($h);
$mpdl->wfits($outfile);

exit;


