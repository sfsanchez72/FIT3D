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


if ($#ARGV<8) {
    print "USE: plot_radial_pt_rss.pl PT.txt INDEX.out X_C Y_C SCALE_PIXEL N.COLUMN2(0...N) LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [error]\n";
    exit;
}

$infile=$ARGV[0];
$infile1=$ARGV[1];
$x0=$ARGV[2];
$x1=$ARGV[3];
$nc1=$ARGV[4];
$nc2=$ARGV[5];
$label1=$ARGV[6];
$label2=$ARGV[7];
$dev=$ARGV[8];
$n=0;
open(FH,"<$infile");
$line=<FH>;
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$x[$n]=sqrt(($data[1]-$x0)**2+($data[2]-$x1)**2);	
	$n++;
    }
}
close(FH);
#print "center = $x0,$x1\n";
$n=0;
open(FH,"<$infile1");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
#	$x[$n]=$data[$nc1];	
	if ($data[$nc2]!=0) {
	    $y[$n]=$data[$nc1]/$data[$nc2];
	    $e_y[$n]=$data[$nc1+1]/$data[$nc2];
#	$delta[$n]=$y[$n]-$x[$n];
	    $n++;
	}
    }
}
$n=$n-1;
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
if ($#ARGV==10) {
    $x_min=$ARGV[9];
    $x_max=$ARGV[10];
}

if ($#ARGV==12) {
    $x_min=$ARGV[9];
    $x_max=$ARGV[10];
    $y_min=$ARGV[11];
    $y_max=$ARGV[12];
}

$fix=1;
if ($#ARGV==13) {
    $x_min=$ARGV[9];
    $x_max=$ARGV[10];
    $y_min=$ARGV[11];
    $y_max=$ARGV[12];
    $fix=$ARGV[13];
}


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
#pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
pgenv($x_min,$x_max,$y_min,$y_max,0,0);
pglabel("$label1","$label2","");
pgsch(1.7);
#print "$n\n @x\n @y\n";
for ($i=0;$i<$n;$i++) {
    $X=$x[$i];
    $Y=$y[$i];
    $eY=$e_y[$i];


 if (($X>$x_min)&&($X<$x_max)&&($Y>$y_min)&&($Y<$y_max)) {
	#pgerrb(1,1,[$X],[$Y],[$eX],0.2);
	#pgerrb(3,1,[$X],[$Y],[$eX],0.2);
     if ($fix==1) {
	 pgerrb(2,1,[$X],[$Y],[$eY],0.2);
	 pgerrb(4,1,[$X],[$Y],[$eY],0.2);
     }
     pgsci(4);
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



exit;





