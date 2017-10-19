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

if ($#ARGV<2) {
    print "USE: spec_comp.pl INPUT_FILE1 INPUT_FILE2 DEVICE [MIN MAX] [WMIN WMAX] [SHIFT]\n";
    exit;
}

$input=$ARGV[0];
$input2=$ARGV[1];
$dev=$ARGV[2];
$def=0;
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}




$y_min=1e12;
$y_max=-1e12;
$n=0;
open(FH,"<$input");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
	$flux[$n]=$data[2];
	$wave[$n]=$data[1];
    } else {
	$flux[$n]=$data[1];
	$wave[$n]=$data[0];
    }
    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n];
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n];
    }
#    print "$n $wave[$n] $flux[$n]\n";
    $n++;
}
close(FH);
$n2=0;
open(FH,"<$input2");
open(OUT,">spec_comp.out");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
	$flux2[$n2]=$data[2];
	$wave2[$n2]=$data[1];
    } else {
	$flux2[$n2]=$data[1];
	$wave2[$n2]=$data[0];
    }
    if ($flux2[$n2]>$y_max) {
	$y_max=$flux2[$n2];
    }
    if ($flux2[$n2]<$y_min) {
	$y_min=$flux2[$n2];
    }
    if ($flux[$n2]!=0) {
	$norm[$n2]=$flux2[$n2]/$flux[$n2];
    } else {
	$norm[$n2]=1;
    }
    print OUT "$data[0] $wave2[$n2] $norm[$n2]\n";
    $n2++;
}
close(OUT);
close(FH);



$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==6) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $def=1;
}
$shift=0;
if ($#ARGV==7) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $shift=$ARGV[7];
    $def=1;
}

for ($j=0;$j<$n;$j++) {
    $wave[$j]=$wave[$j]+$shift;
}

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}
$nmin=$n;
if ($n2<$n) {
    $nmin=$n2;
}
pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","$input vs $input2");
pgline($n,\@wave,\@flux);
pgsci(2);
pgline($n2,\@wave2,\@flux2);
pgsci(3);
pgsls(2);
pgline($n2,\@wave2,\@norm);
pgclose;
pgend;

exit;
