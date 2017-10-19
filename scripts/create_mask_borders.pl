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

if ($#ARGV<7) {
    print "USE: create_mask_borders.pl position_table.txt CRVAL1 CRVAL2 CDELT1 CDELT2 nx ny mask.fits\n";
    exit;
}

$pos_table=$ARGV[0];
$crval1=$ARGV[1];
$crval2=$ARGV[2];
$cdelt1=$ARGV[3];
$cdelt2=$ARGV[4];
$nx=$ARGV[5];
$ny=$ARGV[6];
$mask_file=$ARGV[7];

$mask=zeroes($nx,$ny);
#$crval1=-190.822692871094;
#$crval2=-144.166793823242;

open(FH,"<$pos_table");
while($line=<FH>) {
    chop($line);    
    @data=split(" ",$line);
    if ($data[0] !~ "C") {
	$x=$data[1];
	$y=$data[2];
	$delta=1.8;
	$i_min=int(($x-$delta-$crval1)/$cdelt1);
	$i_max=int(($x+$delta-$crval1)/$cdelt1);
	$j_min=int(($y-$delta-$crval2)/$cdelt2);
	$j_max=int(($y+$delta-$crval2)/$cdelt2);
	#print "$x $y [$i_min:$i_max,$j_min:$j_max]\n";
	if ($i_min<0) {
	    $i_min=0;
	}
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($i_max>($nx)) {
	    $i_max=$nx;
	}
	if ($j_max>($ny)) {
	    $j_max=$ny;
	}

	for ($i=$i_min;$i<$i_max;$i++) {
	    for ($j=$j_min;$j<$j_max;$j++) {
	#	print "$i,$j =";
		
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
