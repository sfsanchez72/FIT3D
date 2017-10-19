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


if ($#ARGV<13) {
    print "USE: disp_lat.pl ASCII_TABLE.txt N.COLUMN1(0...N)  N.COLUMN2(0...N) LABEL1 LABEL2 DEV MIN MAX YMIN YMAX FIX NPOLY CRVAL CDELT\n";
    exit;
}

$infile=$ARGV[0];
$nc1=$ARGV[1];
$nc2=$ARGV[2];
$label1=$ARGV[3];
$label2=$ARGV[4];
$dev=$ARGV[5];
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if (($line !~ "#")&&($line !~ "LINE")) {
	@data=split(" ",$line);
	$x[$n]=$data[$nc1];	
	$y[$n]=$data[$nc2]-$x[$n];
	$delta[$n]=$y[$n];#-$x[$n];
	if ($line !~ "LINE") {
	    $n++;
	}
    }
}
close(FH);
$slice=pdl(@x);
($mean,$rms,$median,$min,$max) = stats($slice);
$x_min=$min-0.2*($median-$min);
$x_max=$max+0.2*($max-$median);
if ($median>2000000) {
    $label1=$label1."-".$median;
    for ($i=0;$i<$n;$i++) {
	$x[$i]=$x[$i]-$median;
    }
    $x_min=$x_min-$median;
    $x_max=$x_max-$median;
}
$slice=pdl(@y);
($mean,$rms,$median,$min,$max) = stats($slice);
$y_min=0.99*$min;
$y_max=1.01*$max;
$NN=$n/int($nbin/6+1);


$slice=pdl(@delta);
($mean,$rms,$median,$min,$max) = stats($slice);
#print "@delta\n";
print "Y-X = $mean $rms $median $min $max\n";
if ($#ARGV==7) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
}

if ($#ARGV==9) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
}

$fix=0;
if ($#ARGV==10) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
}

$npoly=6;
$crval=3745;
$cdelt=2;
if ($#ARGV==13) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $npoly=$ARGV[11];
    $crval=$ARGV[12];
    $cdelt=$ARGV[13];
}


$N=0;
my @xx;
my @yy;
for ($i=0;$i<$n;$i++) {
    if (abs($y[$i]-$mean)<2*$rms) {
	$xx[$NN]=$x[$i];
	$yy[$NN]=$y[$i];
#	print "$xx[$NN] $yy[$NN]  $NN\n";
	if (($yy[$NN]>$y_min)&&($yy[$NN]<$y_max)) {
	    $NN++;
	}
    }
}
$yy[$NN-1]=0;
#$yy[$NN-2]=0;

#print "NN=$NN\n";
#if ($NN>100) {
#print "$npoly $NN\n";

($s_x,$coeff) = fitpoly1d(pdl(@xx),pdl(@yy),$npoly);
@c=list($coeff);
open(OUT,">lat.txt");
print OUT "# Look Up Table\n";
for ($i=0;$i<1970;$i++) {
    $wave[$i]=$crval+$cdelt*$i;
    $out_spec[$i]=0;
    if ($i<11777) {
	for ($k=0;$k<$npoly;$k++) {
	    $out_spec[$i]=$out_spec[$i]+$c[$k]*(($wave[$i])**$k);
	}
    } else {
	$out_spec[$i]=$out_spec[$i-1];
    }
    $out_pix[$i]=$out_spec[$i]/$cdelt;
#    print OUT "$i $wave[$i] $out_spec[$i]\n";
    print OUT "$i $wave[$i] $out_pix[$i]\n";
}
close(OUT);
#}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
#$rms_now=apr($rms);
$rms_now=$rms;
pglabel("$label1","$label2 rms~$rms_now","");
pgsch(2.5);
#print "$n\n @x\n @y\n";
open(OUT,">table_delta_plot.list");
for ($i=0;$i<$NN;$i++) {
    $X=$xx[$i];
    $Y=$yy[$i];
    $D[$i]=$Y;
    print OUT "$i $X $Y\n";
    pgsci(4);
    pgpoint(1,[$X],[$Y],17);
    pgsci(1);
    pgpoint(1,[$X],[$Y],22);
}
#pgpoint($n,\@x,\@y,16);
pgsci(2);
pgline(1970,\@wave,\@out_spec);
pgclose;
pgend;
close(OUT);

$med=median(@D);
$sig=sigma(@D);
$nn=0;
for ($i=0;$i<$n;$i++) {
    if (abs($D[$i]-$med)<2*$sig) {
	$DD[$nn]=$D[$i];
	$nn++;
    }
}
$med=median(@DD);
$sig=sigma(@DD);
print "DELTA X-Y =$med +- $sig\n";

open(OUT,">table_delta_plot.out");
print OUT "$med $sig\n";
close(OUT);

exit;





