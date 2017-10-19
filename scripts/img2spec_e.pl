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

$R3DPATH=$ENV{"R3DPATH"};
$mypath=$R3DPATH."/my.pl"; 
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "$mypath";


$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: img2spec.pl INPUT_FILE.FITS NY OUTPUTFILE.txt\n";
    exit;
}

$input=$ARGV[0];
$e_input="e_".$ARGV[0];
$NY=$ARGV[1];
$output=$ARGV[2];


$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
$pdl_e=rfits("$e_input");
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};
if ($cdelt==0) {
    $cdelt=1;
}
open(FH,">$output");
for ($i=0;$i<$nx;$i++) {
    $k=$i+1;
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($i,$NY);
    if ($flux[$i] eq "BAD") {
	$flux[$i]=0;
    }
    $e_flux[$i]=$pdl_e->at($i,$NY);
    if ($e_flux[$i] eq "BAD") {
	$e_flux[$i]=1e12;
    }
       if ($e_flux[$i]>sqrt(abs($flux[$i]))) {
           $e_flux[$i]=sqrt(abs($flux[$i]));
       }


    print FH "$k $wave[$i] $flux[$i] $e_flux[$i]\n";
}
close(FH);

exit;
