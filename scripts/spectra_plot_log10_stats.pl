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
    print "USE: spectra_plot_log10_stats.pl SET_INPUTFILES DEVICE [MIN MAX] [WMIN WMAX] [YLABEL]\n";
    exit;
}

$sinput=$ARGV[0];
$dev=$ARGV[1];
$def=0;
if ($#ARGV==3) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $def=1;
}




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
print "$nf files found\n";
$wmin=1e12;
$wmax=-1e12;
for ($i=0;$i<$nf;$i++) {
    $input=$file[$i];
#    print "$input ";
    $n=0;
    open(FH,"<$input");
    while ($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($#data==1) {
	    $val=$data[1]*$ratio[$i];
	    if ($val>0) {
		$flux[$i][$n]=log10($val);
	    } else {
		$flux[$i][$n]=-1000;
	    }
	    $wave[$i][$n]=$data[0];
	    if ($wmin>$data[0]) {
		$wmin=$data[0];
	    }
	    if ($wmax<$data[0]) {
		$wmax=$data[0];
	    }
	} else {
	    $val=$data[2]*$ratio[$i];
	    if ($val>0) {
		$flux[$i][$n]=log10($val);
	    } else {
		$flux[$i][$n]=-1000;
	    }
#	    $flux[$i][$n]=$data[2]*$ratio[$i];	
	    $wave[$i][$n]=$data[1];
	    if ($wmin>$data[1]) {
		$wmin=$data[1];
	    }
	    if ($wmax<$data[1]) {
		$wmax=$data[1];
	    }
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
    $nline[$i]=$n;
    close(FH);
#    print "$n\n";
}
#print "$n\n";
#$wmin=$wave[0][0];
#$wmax=$wave[0][$n-1];

#print "[$wmin,$wmax]\n";
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
    $label=$ARGV[6];
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
pgenv($wmin,$wmax,log10($y_min),log10($y_max),0,20);
#pglabel("Wavelength","Flux","");
pglabel("Wavelength (\\A)",$label,"");

pgsch(0.1);
$c=1;
 $n_med=int($n/2);
for ($i=0;$i<$nf;$i++) {    
    $n=$nline[$i];#=$n;
    print "$file[$i] $n $wave[$i][0] $wave[$i][$n-1]\n";
    for ($j=0;$j<$n;$j++) {
	$wave_now[$j]=$wave[$i][$j];
	$flux_now[$j]=$flux[$i][$j];
	$flux_sum[$j]=$flux_sum[$j]+$flux_now[$j];
	
	if ($i==0) {
	    $flux_min[$j]=$flux_now[$j];
	} else {
	    if ($flux_min[$j]>$flux_now[$j]) {
		$flux_min[$j]=$flux_now[$j];
	    }
	}
    }
    print "$flux_now[int($n/2)] $flux_sum[int($n/2)]\n";
 #   <stdin>;
#    pgsci($c);
    pgsci(1);
    pgsls(2);
    pgline($n,\@wave_now,\@flux_now);
#    pgpoint($n,\@wave_now,\@flux_now,1);
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
    pgsch(1.8);


open(OUT,">ratio.mean");
open(MIN,">ratio.min");
for ($j=0;$j<$n;$j++) {
    $flux_sum[$j]=$flux_sum[$j]/$nf;
    $flux_sum_10=10**($flux_sum[$j]);
    $flux_min_10=10**($flux_min[$j]);
    print OUT "$j $wave_now[$j] $flux_sum_10\n";
    print MIN "$j $wave_now[$j] $flux_min_10\n";
}
close(MIN);
close(OUT);
pgsci(2);
pgsls(1);
pgline($n,\@wave_now,\@flux_sum);

pgsci(1);
pgptxt($wmin*1.05,0.92*log10($y_max),0,0,"");

pgclose;
pgend;

exit;
