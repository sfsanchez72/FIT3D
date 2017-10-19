#!/usr/bin/perl
#
#
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

if ($#ARGV<5) {
    print "USE: create_mask_extrapolated.pl position_table.txt NX NY CRVAL1 CRVAL2 mask.fits [DELTA=1.8]\n";
    exit;
}

$pos_table=$ARGV[0];
$NX=$ARGV[1];
$NY=$ARGV[2];
$crval1=$ARGV[3];
$crval2=$ARGV[4];
$mask_file=$ARGV[5];
$delta=1.8;
if ($#ARGV==6) {
    $delta=$ARGV[6];
}


$mask=zeroes($NX,$NY);
#$crval1=-190.822692871094;
#$crval2=-144.166793823242;

open(FH,"<$pos_table");
while($line=<FH>) {
    chop($line);    
    @data=split(" ",$line);
    if ($data[0] !~ "C") {
	$x=$data[1];
	$y=$data[2];
	$i_min=int($x-$delta-$crval1);
	$i_max=int($x+$delta-$crval1);
	$j_min=int($y-$delta-$crval2);
	$j_max=int($y+$delta-$crval2);
	#print "$x $y [$i_min:$i_max,$j_min:$j_max]\n";
	if ($i_min<0) {
	    $i_min=0;
	}
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($i_max>($NX)) {
	    $i_max=$NX;
	}
	if ($j_max>($NY)) {
	    $j_max=$NY;
	}

	for ($i=$i_min;$i<$i_max;$i++) {
	    for ($j=$j_min;$j<$j_max;$j++) {
		#print "$i,$j =";
		set($mask,$i,$j,1.0);
		$val=$mask->at($i,$j);
		#print "$val\n";
	    }
	}
    }
}
close(FH);

$mask->wfits($mask_file);
exit;
