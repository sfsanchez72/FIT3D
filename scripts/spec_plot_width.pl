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
    print "USE: spec_plot_width.pl INPUT_FILE DELTA_A DEVICE [MIN MAX] [WMIN WMAX]\n";
    exit;
}

$input=$ARGV[0];
$width=$ARGV[1];
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
#    print "$#data\n";
    if ($#data>1) {
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
    $n++;
}
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




$ww=$wave[1]-$wave[0];
if ($width>(2*$ww)) {
    $width_pix=$width/($ww);

    $sigma=sqrt(($width*3/2.345)**2-3.3**2);
    $rsigma=$sigma/$ww;
    @flux_c=@flux;
    $box=int(3*$rsigma);
    if ($box>3) {
	for ($i=$box;$i<($n-$box);$i++) {
	    $norm=0;
	    $flux_c[$i]=0;
	    for ($j=$i-$box;$j<($i+$box);$j++) {
		$gaus=exp(-0.5*((($i-$j)/$rsigma)**2));
		$flux_c[$i]=$flux_c[$i]+$flux[$j]*$gaus;
		$norm=$norm+$gaus;
	    }
	    $flux_c[$i]=$flux_c[$i]/$norm;
	}
    }
    






    $nn=int($n/$width_pix);
    print "$nn\n";
    $kk=0;
    for ($i=0;$i<$n;$i=$i+$width_pix) {
	my @F;
	my @W;
	$k=0;
	for ($j=$i;$j<($i+$width_pix);$j++) {
	    $W[$k]=$wave[$j];
	    $F[$k]=$flux_c[$j];
	    $k++;
	}
	$wave_w[$kk]=mean(@W);
	$flux_w[$kk]=mean(@F);
#	print "$kk $wave_w[$kk] $flux_w[$kk]\n";
	$kk++;
    }


#    @flux_w=mean_box($width_pix,\@flux);
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
#
if ($width>(2*$ww)) {
    #pgsci(2);
    pgline($nn,\@wave_w,\@flux_w);
} else {
    pgline($n,\@wave,\@flux);
}
pgsci(2);
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
