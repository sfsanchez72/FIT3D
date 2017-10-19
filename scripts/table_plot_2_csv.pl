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


if ($#ARGV<7) {
    print "USE: table_plot.pl ASCII_TABLE.txt N.C_X1(0...N)  N.C_Y1(0...N) N.C_X2 N.C_Y2 LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [FIX] [DELTA]\n";
    exit;
}

$infile=$ARGV[0];
$ncx1=$ARGV[1];
$ncy1=$ARGV[2];
$ncx2=$ARGV[3];
$ncy2=$ARGV[4];
$label1=$ARGV[5];
$label2=$ARGV[6];
$dev=$ARGV[7];
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    $line =~ s/ //g;
    if ($line !~ "#") {
	@data=split(",",$line);
	$x[$n]=$data[$ncx1];	
	$y[$n]=$data[$ncy1];
	$x2[$n]=$data[$ncx2];	
	$y2[$n]=$data[$ncy2];
	$delta[$n]=$y[$n]-$x[$n];
	$n++;
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

#print "Y-X = $mean $rms $median $min $max\n";
if ($#ARGV==9) {
    $x_min=$ARGV[8];
    $x_max=$ARGV[9];
}

if ($#ARGV==11) {
    $x_min=$ARGV[8];
    $x_max=$ARGV[9];
    $y_min=$ARGV[10];
    $y_max=$ARGV[11];
}

$fix=0;
if ($#ARGV==12) {
    $x_min=$ARGV[8];
    $x_max=$ARGV[9];
    $y_min=$ARGV[10];
    $y_max=$ARGV[11];
    $fix=$ARGV[12];
}

$Delta=0;
if ($#ARGV==13) {
    $x_min=$ARGV[8];
    $x_max=$ARGV[9];
    $y_min=$ARGV[10];
    $y_max=$ARGV[11];
    $fix=$ARGV[12];
    $Delta=$ARGV[13];
}


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.5);           # Set character height
pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
pglabel("$label1","$label2","");
pgsch(2.5);
#print "$n\n @x\n @y\n";
if ($Delta==0) {
    $lim=1e12;
} else {
    $lim=100*abs($Delta);
}
for ($i=0;$i<$n;$i++) {
    $X=$x[$i];
    $Y=$y[$i]-$Delta;
    $D[$i]=$X-$Y;
    if ($D[$i]>1.5) {
	$Y=$y[$i]+2*$Delta;
	$D[$i]=$X-$Y;
    }
#    print "$lim $X $Y $Delta\n";
    if (abs($D[$i])<$lim) {
	pgsci(4);
	pgpoint(1,[$X],[$Y],17);
	pgsci(1);
	pgpoint(1,[$X],[$Y],22);
	$X=$x2[$i];
	$Y=$y2[$i];
	pgsch(2.5);
	pgsci(2);
	pgpoint(1,[$X],[$Y],16);
	pgsci(1);
	pgsch(2.1);
	pgpoint(1,[$X],[$Y],6);
	$KKK++;
    }
}
pgsch(2.5);
#pgpoint($n,\@x,\@y,16);
pgclose;
pgend;


$med=median(@D);
$sig=sigma(@D);
print "DELTA X-Y =$med +- $sig $KKK/$n\n";
exit;





