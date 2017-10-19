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

$galfit="/work1/ssa/galfit/galfit";
$cl="/work1/ssa/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: plot_out_fit.pl INPUT_FILE DEVICE [MIN MAX] [WMIN WMAX]\n";
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
    $wave[$n]=$data[0];
    $flux[$n]=$data[1];
    $flux_mod[$n]=$data[2];
    $flux_res[$n]=$data[3];
    $flux_cont[$n]=$data[4];
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
pgsls(1);
pgsci(3);
pgline($n,\@wave,\@flux_res);
pgsci(1);
pgline($n,\@wave,\@flux);
pgsci(2);
pgline($n,\@wave,\@flux_mod);
pgsls(1);
pgsci(5);
pgline($n,\@wave,\@flux_cont);
pgsls(1);
pgclose;
pgend;

exit;
