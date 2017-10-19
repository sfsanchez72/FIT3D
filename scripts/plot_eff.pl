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

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

$gain=1;
$diam=1;
$factor=1;
if ($#ARGV<1) {
    print "USE: plot_eff.pl RATIO DEVICE [MIN MAX] [WMIN WMAX] [GAIN] [Diameter_tel in m] [Correction_factor]\n";
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
$n=0;
open(FH,"<$input");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
#    print "$#data\n";
    if ($#data>1) {
	$flux[$n]=$data[2];
	$wave[$n]=$data[1];
    } else {
	$flux[$n]=$data[1];
	$wave[$n]=$data[0];
    }
#    print "$wave[$n] $flux[$n]\n";

    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n];
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n];
    }
    $n++;
}
close(FH);
$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==5) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $def=1;
}

if ($#ARGV==6) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $gain=$ARGV[6];
    $def=1;
}

if ($#ARGV==7) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $gain=$ARGV[6];
    $diam=$ARGV[7];
    $def=1;
}

if ($#ARGV==8) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $gain=$ARGV[6];
    $diam=$ARGV[7];
    $factor=$ARGV[8];
    $def=1;
}


$area=3.1415*($diam*0.5*100)**2;


#print "$area\n";

$cdelt=$wave[1]-$wave[0];
for ($i=0;$i<$n;$i++) {
    if ($flux[$i]>0) {
	$eff[$i]=$gain*(1.98644*1e8/$wave[$i])/($flux[$i]*$area)/$factor/$cdelt;
    } else {
	$eff[$i]=0;
    }
}

$n_med=int($n/2);
$y_min=$gain*(1.98644*1e8/$wave[$n_med])/($y_min*$area)/$factor/$cdelt;
$y_max=$gain*(1.98644*1e8/$wave[$n_med])/($y_max*$area)/$factor/$cdelt;

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength (\\A)","Instrument Efficiency","");
pgslw(3);
pgsci(15);
pgline($n,\@wave,\@eff);
pgclose;
pgend;

my @sec;
$nn=0;
for ($i=0;$i<$n;$i++) {
    if (($wave[$i]>$wmin)&&($wave[$i]<$wmax)) {
	$sec[$nn]=$eff[$i];
	$nn++;
    }
}

$median=median(@sec);
$mean=mean(@sec);
$sig=sigma(@sec);

$flux=$mean*$n*($wave[1]-$wave[0]);

print "$median +- $sig (F=$flux)\n";

exit;
