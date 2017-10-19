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


if ($#ARGV<5) {
    print "USE: table_plot.pl ASCII_TABLE.txt N.COLUMN1(0...N)  N.COLUMN2(0...N) LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [FIX] [ONE-to-ONE_LINE 0/1] [SYM] [TOP]\n";
    exit;
}

$infile=$ARGV[0];
$nc1=$ARGV[1];
$nc2=$ARGV[2];
$label1=$ARGV[3];
$label2=$ARGV[4];
$dev=$ARGV[5];

$sym=0;
$fix=0;
$one=0;
$top=0;
if ($#ARGV==12) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
    $sym=$ARGV[12];
}



if ($#ARGV==13) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
    $sym=$ARGV[12];
    $top=$ARGV[13];
}

$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(",",$line);
	$x[$n]=$data[$nc1];	
	$y[$n]=$data[$nc2];
	$ex[$n]=0.5*$data[$nc1+1];	
	$ey[$n]=0.5*$data[$nc2+1];
	$delta[$n]=$y[$n]-$x[$n];
	$SYM[$n]=$data[$sym];
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


if ($#ARGV==10) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
}


if ($#ARGV==11) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
}

if ($#ARGV==12) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
    $sym=$ARGV[12];
}

if ($#ARGV==13) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
    $sym=$ARGV[12];
    $top=$ARGV[13];
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
if ($top==0) {
    pglabel("$label1","$label2","");
} else {
    pglabel("","$label2","$label1");
}
if ($one==1) {
    pgsci(2);
    pgsls(2);
    pgslw(3);
    pgline(2,[$x_min,$x_max],[$x_min,$x_max]);
    pgsci(1);
    pgsls(1);
    pgslw(1);
}
pgsch(2.5);
#print "$n\n @x\n @y\n";
for ($i=0;$i<$n;$i++) {
    $X=$x[$i];
    $Y=$y[$i];
    $eX=$ex[$i];
    $eY=$ey[$i];

    $D[$i]=$X-$Y;

#    print "$i $X $Y\n";
    if (($X>$x_min)&&($X<$x_max)&&($Y>$y_min)&&($Y<$y_max)) {
	#pgerrb(1,1,[$X],[$Y],[$eX],0.2);
	#pgerrb(3,1,[$X],[$Y],[$eX],0.2);
	pgerrb(2,1,[$X],[$Y],[$eY],0.2);
	pgerrb(4,1,[$X],[$Y],[$eY],0.2);
	pgsci(4+4*$SYM[$i]);
	if ($SYM[$i]==0) {
	    pgpoint(1,[$X],[$Y],17);
	    pgsci(1);
	    pgpoint(1,[$X],[$Y],22);
	} else {
	    pgsch(2.5);
	    pgpoint(1,[$X],[$Y],16);
	    pgsci(1);
	    pgsch(2);
	    pgpoint(1,[$X],[$Y],6);
	    pgsch(2.5);
	}
    }
}
#pgpoint($n,\@x,\@y,16);
pgclose;
pgend;


$med=median(@D);
$sig=sigma(@D);
print "DELTA X-Y =$med +- $sig (N=$n)\n";
#$med=median(@R);
#$sig=sigma(@R);
#print "RAT X/Y =$med +- $sig (N=$n)\n";

#$slice=pdl(@delta);
#($mean,$rms,$median,$min,$max) = stats($slice);
#print "($mean,$rms,$median,$min,$max)\n";

exit;





