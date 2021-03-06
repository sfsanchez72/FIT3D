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


if ($#ARGV<3) {
    print "USE: imarith.pl INPUT1.FITS OPERATOR(+,-,/,*) INPUT_SPEC.TXT OUTPUT.FITS\n";
    exit;
}

$infile1=$ARGV[0];
$operator=$ARGV[1];
$infile2=$ARGV[2];
$outfile=$ARGV[3];

$a_in1=rfits($infile1);

$h=$a_in1->gethdr;
($nx,$ny)=$a_in1->dims;


$n=0;
open(FH,"<$infile2");
while($line=<FH>) {
    @data=split(" ",$line);
    $index[$n]=$data[0];
    $wave[$n]=$data[1];
    $flux_ini[$n]=$data[2];
    $n++;
}
close(FH);

if ($nx!=$n) {
    print "Image and spectra size do not match ($nx!=$n)!\n";
    exit;
}

$a_in2=zeroes($nx);

for ($i=0;$i<$nx;$i++) {
    set($a_in2,$i,$flux_ini[$i]);
}

for ($j=0;$j<$ny;$j++) {
    $jend=$j+1;
    $t=$a_in1->slice(":,$j");

    if ($operator eq "+") {
	$t .=$t+$a_in2;
    }

    if ($operator eq "-") {
	$t .=$t-$a_in2;
    }

    if ($operator eq "/") {
	$t .=$t/$a_in2;
    }

    if ($operator eq "*") {
	$t .=$t*$a_in2;
    }
    
}

$a_in1->sethdr($h);
$a_in1->wfits($outfile);

exit;
