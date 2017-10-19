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

if ($#ARGV<1) {
    print "USE: spec_plot.pl INPUT_FILE DEVICE  [MIN MAX] [WMIN WMAX]\n";
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
    if ($line !~ "#") {
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

$crval=$wave[0];
$cdelt=($wave[1]-$wave[0]);
$npix=int(($wave[$n-1]-$crval)/$cdelt);
for ($j=0;$j<$npix;$j++) {
    $wave_new[$j]=$crval+$cdelt*$j;
}

$out_spec_pdl = interpol(pdl(@wave_new), pdl(@wave), pdl(@flux));
open(OUT,">spec_resample.out");
for ($j=0;$j<$npix;$j++) {
    $val=$out_spec_pdl->at($j);
    print OUT "$j $wave_new[$j] $val\n";
}
close(OUT);


exit;
