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

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<6) {
    print "USE: DAR_estimate.pl X_Y_mesurements.txt NOUT X Y NP_X NP_Y DAR_estimation.txt\n";
    exit;
}

$input=$ARGV[0];
$nz=$ARGV[1];
$x0=$ARGV[2];
$y0=$ARGV[3];
$npoly_x=$ARGV[4];
$npoly_y=$ARGV[5];
$output=$ARGV[6];


$NZ=0;
open(FH,"<$input");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $point[$NZ]=$data[0];
    $x[$NZ]=$data[1]-$x0;
    $y[$NZ]=$data[2]-$y0;
    $NZ++;
}



$point_pdl=pdl(@point);
$y_pdl = pdl(@y);
($s_y,$coeff) = fitpoly1d $point_pdl,$y_pdl,$npoly_y;

for ($j=0;$j<$nz;$j++) {
    $y_out[$j]=0;
    for ($i=0;$i<$npoly_y;$i++) {
	$C=$coeff->at($i);
	$y_out[$j]=$y_out[$j]+$C*($j**$i);
    }
}

$x_pdl = pdl(@x);
($s_x,$coeff) = fitpoly1d $point_pdl,$x_pdl,$npoly_x;

for ($j=0;$j<$nz;$j++) {
    $x_out[$j]=0;
    $id[$j]=$j;
    for ($i=0;$i<$npoly_x;$i++) {
	$C=$coeff->at($i);
	$x_out[$j]=$x_out[$j]+$C*($j**$i);
    }
}

open(OUT,">$output");
for ($i=0;$i<$nz;$i++) {
    print OUT "$i $x_out[$i] $y_out[$i]\n";
}
close(OUT);

($xmin,$xmax)=minmax(@x_out);
($ymin,$ymax)=minmax(@y_out);

if ($xmin<$ymin) {
    $min=$xmin;
} else {
    $min=$ymin;
}

if ($xmax>$ymax) {
    $max=$xmax;
} else {
    $max=$ymax;
}

pgbegin(0,"/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
pgenv(1,$nz,$min-1,$max+1,0,0);
pgsci(2);
pgpoint($NZ,\@point,\@x,3);
pgsci(4);
pgpoint($NZ,\@point,\@y,3);
pgsci(1);
pgline($nz,\@id,\@x_out);
pgline($nz,\@id,\@y_out);
pgclose;
pgend;

exit;
