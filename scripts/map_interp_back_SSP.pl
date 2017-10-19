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




if ($#ARGV<5) {
    print "USE: map_interp_back.pl slice.txt/POS_TABLE out.model PREFIX_OUT INTERP_MODE GRID_FUNC  GRID_PIX \n";
    exit;
}

$slice=$ARGV[0];
$out_model=$ARGV[1];
$prefix=$ARGV[2];
$int_mode=$ARGV[3];
$int_opt=$ARGV[4];
$int_dpix=$ARGV[5];
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
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
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

@NAME=("n","i","j","k","chi","age","met","Av_ssp","redshift","sigma","flux","red_abs","red_em","med_flux","sig_flux","age_mass","met_mass");

$n_row=0;
open(FH,"<$out_model");
while($line=<FH>) {
    chop($line);    
    @data=split(" ",$line);
    $NDATA=$#data+1;
    for ($i=0;$i<$NDATA;$i++) {
	$DATA[$n_row][$i]=$data[$i];	
    }    
    $n_row++;
}
close(FH);

print "$NDATA $n_row\n";
for ($i=4;$i<$NDATA;$i++) {
    $outfile=$prefix."_".$NAME[$i].".fits";
    $slice=$prefix."_".$NAME[$i].".slice";
    open(JUNK,">junk_int.tmp");
    for ($j=0;$j<$n_row;$j++) {
	printf JUNK "$j $x[$j] $y[$j] $DATA[$j][$i]\n";
    }
    close(JUNK);
    system("rm $outfile_f");
    $call="interpol -if junk_int.tmp -of ".$outfile." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
    system($call);
    $call="cp junk_int.tmp ".$slice;
    system($call);
}

exit;
