#!/usr/bin/perl

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

if ($#ARGV<3) {
    print "USE: create_mask_list.pl NX NY list.txt map.fits\n";
    print "list.txt\n";
    print "X Y radius\n";
    exit;
}

$nx=$ARGV[0];
$ny=$ARGV[1];
$input_file=$ARGV[2];
$output_file=$ARGV[3];
$pdl_output=ones($nx,$ny);
$n=0;
open(FH,"<$input_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $ic=$data[0];
    $jc=$data[1];
    $rc=$data[2];
    $i1=$ic-1.5*$rc;
    $i2=$ic+1.5*$rc;
    $j1=$jc-1.5*$rc;
    $j2=$jc+1.5*$rc;
    if ($i1<0) {
	$i1=0;
    }
    if ($i2>$nx) {
	$i2=$nx;
    }
    if ($j1<0) {
	$j1=0;
    }
    if ($j2>$ny) {
	$j2=$ny;
    }
    for ($i=$i1;$i<$i2;$i++) {
	for ($j=$j1;$j<$j2;$j++) {
	    $dist=sqrt(($i-$ic)**2+($j-$jc)**2);
	    if ($dist<=$rc) {
		set($pdl_output,$i,$j,0);
	    }
	}
    }

}
close(FH);



$pdl_output->wfits($output_file);
exit;


