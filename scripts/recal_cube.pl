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


if ($#ARGV<2) {
    print "USE: recal_cube.pl CUBE.fits RATIO.spec OUTPUTCUBE.fits\n";
    exit;
}

$infile=$ARGV[0];
$ratio=$ARGV[1];
$outfile=$ARGV[2];


$pdl=rfits($infile);
($nx,$ny,$nz)=$pdl->dims;


$n=0;
open(FH,"<$ratio");
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
	$index[$n]=$data[0];
	$wave[$n]=$data[1];
	$rat[$n]=$data[2];
	$n++;
    }
}
close(FH);

if ($n!=$nz) {
    print "$nz!=$n. Not possible to scale\n";
}


$pdl_rat=pdl(@rat);

for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$t=$pdl->slice("($i),($j),");
	$t .= $t*$pdl_rat;
    }
}

$pdl->wfits($outfile);

