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

if ($#ARGV<3) {
    print "USE: plot_CAFE_order.pl SET_INPUTFILES DEVICE MIN.ORDER MAX.ORDER [MIN MAX] [WMIN WMAX]\n";
    exit;
}

$sinput=$ARGV[0];
$dev=$ARGV[1];
$min_o=$ARGV[2];
$max_o=$ARGV[3];
$def=0;
if ($#ARGV==5) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $def=1;
}
$NY=$max_o-$min_o;



$y_min=1e12;
$y_max=-1e12;
$nf=0;
open(C,"<$sinput");
while ($line=<C>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==0) {
	$file[$nf]=$data[0];
	$ratio[$nf]=1;
    } else {
	$file[$nf]=$data[0];
	$ratio[$nf]=$data[1];
    }
#    print "$ratio[$nf]\n";
    $nf++;
}
close(C);
#print "$nf files found\n";
$wmin=1e12;
$wmax=-1e12;
for ($i=0;$i<$nf;$i++) {
    $input=$file[$i];
    #print "$input ";
    $n=0;
    open(FH,"<$input");
    while ($line=<FH>) {
	chop($line);
	if ($line !~ "#") {
	    @data=split(" ",$line);
	    if ($#data==2) {
		$flux[$i][$n]=$data[2]*$ratio[$i];	
		$wave[$i][$n]=$data[1];
		if ($wmin>$data[1]) {
		    $wmin=$data[1];
	    }
		if ($wmax<$data[1]) {
		    $wmax=$data[1];
		}
	    } else {
		$flux[$i][$n]=$data[1]*$ratio[$i];	
		$wave[$i][$n]=$data[0];
		
		if ($wmin>$data[0]) {
		    $wmin=$data[0];
		}
		if ($wmax<$data[0]) {
		    $wmax=$data[0];
		}
	    }
#	print "$flux[$i][$n] $wave[$i][$n] $ratio[$i]\n";

	    if (($i>=$min_o)&&($i<$max_o)) {
		
		if ($flux[$i][$n]>$y_max) {
		    $y_max=$flux[$i][$n];
		}
		if ($flux[$i][$n]<$y_min) {
		    $y_min=$flux[$i][$n];
		}
	    }
	$n++;
	}
    }
    $nline[$i]=$n;
    close(FH);
#    print "$n\n";
}
#print "$n\n";
#$wmin=$wave[0][0];
#$wmax=$wave[0][$n-1];
#print "[$wmin,$wmax]\n";
if ($#ARGV==7) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $wmin=$ARGV[6];
    $wmax=$ARGV[7];
    $def=1;
}

$y_min=$y_min-0.2*abs($y_min);
#print "$y_min\n";
if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgsch(10/$NY);           # Set character height
pgsubp(1,$NY);
$c=1;
for ($i=$min_o;$i<$max_o;$i++) {    
    my @wave_now;
    my @flux_now;
    
    for ($j=0;$j<$n;$j++) {
	$wave_now[$j]=$wave[$i][$j];
	$flux_now[$j]=$flux[$i][$j];
    }
    $W1=$wave[$i][0];
    $W2=$wave[$i][$n-1];
    $med=median(@flux_now);
    $sig=sigma(@flux_now);
    $y_min=$med-5*$sig;
    $y_max=$med+5*$sig;
#    print "$i $n $W1 $W2\n";
    pgenv($W1,$W2,$y_min,$y_max,0,0);
    $order=$i+59;
    $n=$nline[$i];#=$n;
    #print "$file[$i] $n $wave[$i][0] $wave[$i][$n-1] ";
    $wave_mean=0;
    $sum_flux=0;

    pgsci(2);    
    pgline($n,\@wave_now,\@flux_now);
    $kk=$i+1;    
    srandom;
    pgsch(3);
    pgsci(14);
    pgptxt($W1+0.02*($W2-$W1),$y_max-0.2*($y_max-$y_min),0,0,"$W1-$W2 ORDER $order");
    pgsch(10/$NY);
    pgsci(1);
}
pgclose;
pgend;

exit;
