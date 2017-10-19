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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<4) {
    print "USE: redshift_find.pl elines_find.txt elines_restframe.txt delta_AA min_redshift max_redshift\n";
    exit;
}

$elines_find=$ARGV[0];
$elines_rest=$ARGV[1];
$dA=$ARGV[2];
$zmin=$ARGV[3];
$zmax=$ARGV[4];

$nf=0;
open(FH,"<$elines_find");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$line_f[$nf]=$data[0];
	$nf++;
    }
}
close(FH);

$nr=0;
open(FH,"<$elines_rest");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$line_r[$nr]=$data[0];
	$nr++;
    }
}
close(FH);

$z=$zmin;
$nmatch_max=0;
while($z<$zmax) {
    $nmatch=0;
    my @Zf;
    my @line_rf;
#    for ($i=0;$i<$nr;$i++) {
#	$line_rf[$i]=$line_r[$i]*(1+$z);
#    }

    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$nr;$i++) {
	    $line_rf=$line_r[$i]*(1+$z);
	    if (abs($line_f[$j]-$line_rf)<$dA) {
		$Zf[$nmatch]=$line_f[$j]/$line_r[$i]-1;
		$nmatch++;
	    }
	}
    }
    if ($nmatch>$nmatch_max) {
	$Z=mean(@Zf);
	$eZ=sigma(@Zf);
#	print "$Z $eZ $nmatch\n";
	$nmatch_max=$nmatch;
    }
#    print "* $z\n";
    $z=$z+0.001;
}

print "$Z $eZ $nmatch_max/$nf\n";

open(OUT,">redshift_find.out");
print OUT "$Z $eZ $nmatch_max/$nf\n";
close(OUT);

exit;
