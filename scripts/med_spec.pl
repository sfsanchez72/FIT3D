#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<3) {
    print "USE: med_spec.pl intput_spec.txt output_spec.txt WIDTH DEV\n";
    exit;
}

$input=$ARGV[0];
$outfile=$ARGV[1];
$box=$ARGV[2];
$dev=$ARGV[3];

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
#$med=abs(median(@flux));
#$y_max=3*$med;
#$y_min=-0.5*$med;
#print "$med $y_min $y_max\n";
if ($#ARGV==7) {
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $y_min=$ARGV[6];
    $y_max=$ARGV[7];
}

#print "$#ARGV\n";
if ($#ARGV==5) {
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
}


@med_spec=median_filter($box,\@flux);

open(OUT,">$outfile");
for ($i=0;$i<$n;$i++) {
    $val=$flux[$i]-$med_spec[$i];
    print OUT "$i $wave[$i] $val\n";
}
close(OUT);

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
pgline($n,\@wave,\@flux);
pgsci(8);
pgline($n,\@wave,\@med_spec);
pgclose;
pgend;


exit;


