#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
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


if ($#ARGV<3) {
    print "USE: dust_spec.pl INPUT_FILE Av redshift OUTPUT_FILE [RV=3.1]\n";
    exit;
}

$input=$ARGV[0];
$Av=$ARGV[1];
$redshift=$ARGV[2];
$output=$ARGV[3];
$Rv=3.1;
if ($#ARGV==4) {
    $Rv=$ARGV[4];
}

$n=0;
open(FH,"<$input");
open(OUT,">$output");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data>1) {
	$flux[$n]=$data[2];
	$wave[$n]=$data[1];
    } else {
	$flux[$n]=$data[1];
	$wave[$n]=$data[0];
    }
    $wave_res=$wave[$n]/(1+$redshift);
    $dust_rat=A_l($Rv,$wave_res);
    $dust=10**(-0.4*$Av*$dust_rat);  
    $val=$flux[$n]*$dust;
    print OUT "$n  $wave[$n] $val\n";
    $n++;
}
close(OUT);
close(FH);

open(LIST,">dust_spec.list");
print LIST "$input\n";
print LIST "$output\n";
close(LIST);


exit;
