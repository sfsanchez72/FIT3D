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

if ($#ARGV<6) {
    print "USE: spec_plot.pl SET_INPUTFILES DEVICE MIN MAX WMIN WMAX NAME\n";
    exit;
}

$sinput=$ARGV[0];
$dev=$ARGV[1];
$def=0;
$min=$ARGV[2];
$max=$ARGV[3];
$wmin=$ARGV[4];
$wmax=$ARGV[5];
$name=$ARGV[6];
$def=1;





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
    print "$ratio[$nf]\n";
    $nf++;
}
close(C);
print "$nf files found\n";
for ($i=0;$i<$nf;$i++) {
    $input=$file[$i];
    print "$input ";
    $n=0;
    open(FH,"<$input");
    while ($line=<FH>) {
	chop($line);
	if ($line !~ "#") {
	    @data=split(" ",$line);
	    if ($#data==2) {
		$flux[$i][$n]=$data[2]*$ratio[$i];	
		$wave[$i][$n]=$data[1];
	    } else {
		$flux[$i][$n]=$data[1]*$ratio[$i];	
		$wave[$i][$n]=$data[0];
		
	    }
#	print "$flux[$i][$n] $wave[$i][$n] $ratio[$i]\n";

	if ($flux[$i][$n]>$y_max) {
	    $y_max=$flux[$i][$n];
	}
	if ($flux[$i][$n]<$y_min) {
	    $y_min=$flux[$i][$n];
	}
	$n++;
	}
    }
    $nline[$i]=$n;
    close(FH);
    print "$n\n";
}
print "$n\n";
#$wmin=$wave[0][0];
#$wmax=$wave[0][$n-1];
print "[$wmin,$wmax]\n";

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
pglabel("Wavelength (\\A)","Flux (10\\u-16\\d erg s\\u-1\\d cm\\u-2\\d \\A\\u-1\\d )","");
$c=1;
for ($i=0;$i<$nf;$i++) {    
    $n=$nline[$i];#=$n;
    print "$file[$i] $n $wave[$i][0] $wave[$i][$n-1] ";
    $wave_mean=0;
    $sum_flux=0;
    for ($j=0;$j<$n;$j++) {
	$wave_now[$j]=$wave[$i][$j];
	$flux_now[$j]=$flux[$i][$j];
	if ($flux_now[$j]>0) {
	    $wave_mean=$wave_mean+$flux_now[$j]*$wave_now[$j];
	    $sum_flux=$sum_flux+$flux_now[$j];
	}
    }
    $wave_mean=$wave_mean/$sum_flux;
    print "$wave_mean\n";
    pgsci($c);
    pgline($n,\@wave_now,\@flux_now);
    pgsch(1.8);
    $kk=$i+1;    
    srandom;
    $pos=int($n/5+$i*5);#+random(3);
#    pgsci(1);
    #pgptxt($wave_now[$pos],$flux_now[$pos],0,0,$kk);
    pgsch(1.2);
    $c++;
    if ($c>15) {
	$c=1;
    }
}
pgsci(1);
pgsch(1.4);
pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.1*($y_max-$y_min),0,0,"$name");
pgclose;
pgend;

exit;
