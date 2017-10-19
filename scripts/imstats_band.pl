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


if ($#ARGV<6) {
    print "USE: imstats.pl INPUT.FITS X_C Y_C X_WIDTH Y_WIDTH MIN_VAL MAX_VAL\n";
    exit;
}

$infile=$ARGV[0];
$x_c=$ARGV[1];
$y_c=$ARGV[2];
$dx=$ARGV[3];
$dy=$ARGV[4];
$min=$ARGV[5];
$max=$ARGV[6];
$x_ini=$x_c-$dx;
$x_end=$x_c+$dx;
$y_ini=$y_c-$dy;
$y_end=$y_c+$dy;
$pdl=rfits($infile);
($nx,$ny)=$pdl->dims;
$slice=$pdl->slice("$x_ini:$x_end,$y_ini:$y_end");

($nx1,$ny1)=$slice->dims;
$k=0;
for ($i=0;$i<$nx1;$i++) {
    for ($j=0;$j<$ny1;$j++) {
	$val=$slice->at($i,$j);
	if (($val>$min)&&($val<$max)) {
	    $flux[$k]=$val;
	    $k++;
	}
    }
}
$S=pdl(@flux);

($mean,$rms,$median,$min,$max) = stats($S);

print "$mean $rms $median $min $max\n";

exit;





