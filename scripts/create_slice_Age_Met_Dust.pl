#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";




if ($#ARGV<2) {
    print "USE: create_slice_Age_Met.pl POS_TABLE out_AGE_MET PREFIX_OUT [MIN_FLUX_ALLOWED]\n";
    exit;
}

$slice=$ARGV[0];
$out_model=$ARGV[1];
$prefix=$ARGV[2];
if ($#ARGV==3) {
    $min_limit=$ARGV[3];
} else {
    $min_limit=-1e12;
}
$over=0;

$ns=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;
$dpix=10e10;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=(-1)*$data[1];
	$y[$ns]=(-1)*$data[2];
	if ($x_min>$x[$ns]) {
	    $x_min=$x[$ns];
	}
	if ($y_min>$y[$ns]) {
	    $y_min=$y[$ns];
	}
	if ($x_max<$x[$ns]) {
	    $x_max=$x[$ns];
	}
	if ($y_max<$y[$ns]) {
	    $y_max=$y[$ns];
	}
	if ($ns>0) {
	    if (($dpix>abs($x[$ns]-$x[$ns-1]))&&(abs($x[$ns]-$x[$ns-1])!=0)) {
		$dpix=abs($x[$ns]-$x[$ns-1]);
	    }
	}
	$ns++;
    }
}
close(FH);

$nx=abs($x_max-$x_min)/$dpix+1;
$ny=abs($y_max-$y_min)/$dpix+1;
#print "GRID=$x_min $x_max $y_min $y_max $dpix [$nx,$ny]\n";
#exit;

$n_row=0;
open(FH,"<$out_model");
while($line=<FH>) {
    @data=split(" ",$line);
    $flux[$n_row]=$data[0];
    $val[0][$n_row]=$data[2]; #AGE
    $val[1][$n_row]=$data[3]; #E_AGE
    $val[2][$n_row]=$data[4]; #MET
    $val[3][$n_row]=$data[5]; #E_MET
    $val[4][$n_row]=$data[6]; #DUST
    $val[5][$n_row]=$data[7]; #E_DUST
    $n_row++;
}
close(FH);

$filename[0]=$prefix."_AGE.txt";
$filename[1]=$prefix."_eAGE.txt";
$filename[2]=$prefix."_MET.txt";
$filename[3]=$prefix."_eMET.txt";
$filename[4]=$prefix."_DUST.txt";
$filename[5]=$prefix."_eDUST.txt";
$file_flux=$prefix."_flux.txt";
for ($k=0;$k<6;$k++) {
    open(FILE,">$filename[$k]");
    for ($i=0;$i<$n_row;$i++) {
	if ($flux[$i]>$min_limit) {
	    print FILE "$id[$i] $x[$i] $y[$i] $val[$k][$i]\n";
	} else {
	    print FILE "$id[$i] $x[$i] $y[$i] 0\n";
	}
    }
    close(FILE);
}
open(FILE,">$file_flux");
    for ($i=0;$i<$n_row;$i++) {
	print FILE "$id[$i] $x[$i] $y[$i] $flux[$i]\n";
    }
close(FILE);

exit;
