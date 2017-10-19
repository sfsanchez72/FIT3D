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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<1) {
    print "USE: spec_plot.pl ORDER DEVICE  [MIN MAX] [WMIN WMAX]\n";
    exit;
}

$input=$ARGV[0];
$dev=$ARGV[1];
$def=0;
if ($#ARGV==3) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $def=1;
}




$y_min=1e12;
$y_max=-1e12;
for ($n=0;$n<100;$n++) {
    $wave[$n]=$n/50-1;
    $flux[$n]=legendre($input,$wave[$n]);
    $flux2[$n]=legendre_new($input,$wave[$n]);
#    print "$wave[$n] $flux[$n]\n";
    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n]+0.1;
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n]-0.1;
    }
}


$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==5) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $def=1;
}




if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
pgline($n,\@wave,\@flux);
pgsci(2);
pgpoint($n,\@wave,\@flux2,22);
pgclose;
pgend;

my @sec;
$nn=0;
for ($i=0;$i<$n;$i++) {
    if (($wave[$i]>$wmin)&&($wave[$i]<$wmax)) {
	$sec[$nn]=$flux[$i];
	$nn++;
    }
}

$median=median(@sec);
$mean=mean(@sec);
$sig=sigma(@sec);

$flux=$mean*$n*($wave[1]-$wave[0]);

print "$median +- $sig (F=$flux)\n";

exit;
