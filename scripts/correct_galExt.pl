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

if ($#ARGV<4) {
    print "USE: correct_galExt.pl INPUT_FILE A_lambda wavelength OUTPUT_FILE DEV [MIN MAX] [WMIN WMAX]\n";
    exit;
}

$input=$ARGV[0];
$Al=$ARGV[1];
$lambda=$ARGV[2];
$output=$ARGV[3];
$dev=$ARGV[4];
$def=0;
if ($#ARGV==6) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $def=1;
}

$a_V=A_l(3.1,5505);
$a_l=A_l(3.1,$lambda);

$AV=($a_V/$a_l)*$Al;
#$AV=($a_V/$a_l)*$Al;
#$AV=($a_l)*$Al;


$y_min=1e12;
$y_max=-1e12;
$n=0;
open(FH,"<$input");
open(OUT,">$output");
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
    
    $C=$AV*A_l(3.1,$wave[$n]);
    $flux_cor[$n]=$flux[$n]*(10**(0.4*($C)));
    print OUT "$n $wave[$n] $flux_cor[$n]\n";

    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n];
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n];
    }

    if ($flux_cor[$n]>$y_max) {
	$y_max=$flux_cor[$n];
    }
    if ($flux_cor[$n]<$y_min) {
	$y_min=$flux_cor[$n];
    }


    $n++;
}
    close(OUT);
close(FH);
$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==8) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $wmin=$ARGV[7];
    $wmax=$ARGV[8];
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
pgline($n,\@wave,\@flux_cor);

pgclose;
pgend;

exit;
