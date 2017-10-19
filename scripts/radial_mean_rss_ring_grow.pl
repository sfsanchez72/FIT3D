#!/usr/bin/perl
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
#use PDL::Graphics::TriD;
#use PDL::Graphics::TriD::Image;

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<5) {
    print "USE: radial_sum_rss_ring.pl RSS.fits POS_TABLE.txt Delta_R X_C Y_C OUTPUT.RSS.fits [NMIN]\n";
    exit;
}

$input=$ARGV[0];
$slice=$ARGV[1];
$Dr=$ARGV[2];
$x_c=$ARGV[3];
$y_c=$ARGV[4];
$output=$ARGV[5];
$e_output="e_".$output;
if ($#ARGV==6) {
    $nmin=$ARGV[6];
}

$r_max=0;
$n=0;
open(FH,"<$slice");
$line=<FH>;
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $x[$n]=$data[1];
    $y[$n]=$data[2];
    $flux[$n]=$data[3];
    $r[$n]=sqrt(($x[$n]-$x_c)**2+($y[$n]-$y_c)**2);   
   
    $nk=$r[$n]/$Dr;
    $c[$n]=2+int($nk/2);
    $s[$n]=3+int($nk/2);
    $S[$n]=2.5-$nk/15;

    if ($r_max<$r[$n]) {
	$r_max=$r[$n];	
    }
    $n++;
}
close(FH);

$nr=int($r_max/$Dr)+1;

$pdl_in=rfits("$input");
($nx,$ny)=$pdl_in->dims();
$pdl_out=zeroes($nx,$nr);
$h=$pdl_in->gethdr;

for ($i=0;$i<$nr;$i++) {
    $r_min=$Dr*$i**1.1;
    $r_max=$Dr*($i+1)**1.2;
    $r_mean=($r_min+$r_max)*0.5;
    $nsum_old=$nsum;
    $nsum=0;
    $sum_all=0;
    $t=$pdl_out->slice(":,($i)");
    for ($j=0;$j<$ny;$j++) {
	#print "$j $r[$j]=sqrt(($x[$j]-$x_c)**2+($y[$j]-$y_c)**2);    $r_min $r_max\n";
      if (($r[$j]>$r_min)&&($r[$j]<=$r_max)) {
	    $slice=$pdl_in->slice(":,($j)");
	    $suma=sum($slice);
	    if (($suma!=0)&&($suma ne "nan")) {
	      $t.=$t+$slice;	    
	      $nsum++;
	    }
	  }
    }
    $area[$i]=$nsum;
    if ($nsum==0) {
	$nr=$i;
    }
    print "$i $r_min $r_max $nsum $area[$i] $r_mean\n";
}

for ($i=0;$i<$nr;$i++) {
  $t=$pdl_out->slice(":,($i)");
  $t.=$t/$area[$i];
}


$$h{"NAXIS2"}=$nr;
$pdl_out->sethdr($h);
$pdl_out->wfits($output);

exit;
