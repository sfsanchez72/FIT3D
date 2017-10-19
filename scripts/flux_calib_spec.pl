#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
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
use PDL::Image2D;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<5) {
    print "USE: flux_calib.pl NONCALIBRATED.FITS RATIO.TXT EXPTIME EXTINC_AV AIRMASS OUTPUTFILE.FITS \n";
    exit;
}
$infile=$ARGV[0];
$ratio=$ARGV[1];
$exptime=$ARGV[2];
$ext_AV=$ARGV[3];
$airmass=$ARGV[4];
$outfile=$ARGV[5];

$pdl=rfits($infile);
$h=$pdl->gethdr;
($nx,$ny)=$pdl->dims;


$n=0;
open(FH,"<$ratio");
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
	$index[$n]=$data[0];
	$wave[$n]=$data[1];
	$Krc=0.085*($wave[$n]/5450)**(-4);
#	$ext=1.1*$Krc+(0.8*$ext_AV*$airmass-1.1*0.085)*($wave[$n]/5450)**(-0.8);
	$ext=1.1*$Krc+(0.8*$ext_AV-1.1*0.085)*($wave[$n]/5450)**(-0.8);
	$Fext=10**(0.4*$ext*$airmass);
#	$flux_unc_ini[$n_unc]=$data[2]*$factor/($Fext);
	$DIV=$Fext*$exptime;
	if ($DIV>0) {
	    $rat[$n]=$data[2]/($Fext*$exptime);
	} else {
	    $rat[$n]=$data[2];
	}
	$n++;
    }
}
close(FH);

if ($nx!=$n) {
    print "The number of rows in the RATIO.TXT file does not match\n";
    print " the pixels in the FITSFILE.FITS ($nx!=$n)\n";
#   exit;
    if ($nx<$n) {
	exit;
    } else {
	for($i=0;$i<$n;$i++) {
	    $t = $pdl->slice("$i,:");
	    $t .= $t*$rat[$i];
	}
	$k=$n-1;
	$$h{NAXIS1}=$n;
	$pdl_sec=$pdl->slice("0:$k,:");
	$pdl_sec->sethdr($h);
	$pdl_sec->wfits($outfile);
    }

} else  {
    for($i=0;$i<$nx;$i++) {
	$t = $pdl->slice("$i,:");
	$t .= $t*$rat[$i];
    }
    $pdl->wfits($outfile);
}




exit;
