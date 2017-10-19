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
use PDL::Image2D;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<3) {
    print "USE: mask_spec.pl input_spec.txt mask_list output_spec.txt wave_column\n";
    print "mask_list: Start_wave End_wave [to remove]\n";
    exit;
}

$in_spec=$ARGV[0];
$mask_file=$ARGV[1];
$out_spec=$ARGV[2];
$wave_column=$ARGV[3];

$n_mask=0;
open(FH,"<$mask_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $start_mask[$n_mask]=$data[0];
    $end_mask[$n_mask]=$data[1];
    $n_mask++;
}
close(FH);
print "$n_mask spectral regions to mask\n";


$n=0;
open(FH,"<$in_spec");
while($line=<FH>) {
    @data=split(" ",$line);
    $wave[$n]=$data[$wave_column];
    $flux[$n]=$data[$wave_column+1];
    $n++;
}
close(FH);

for ($i=0;$i<$n;$i++) {
    $copy=1;
    for ($j=0;$j<$n_mask;$j++) {
	if (($wave[$i]>$start_mask[$j])&&($wave[$i]<$end_mask[$j])) {
	    $copy=0;
	}
    }
    if ($copy==1) {
	$wave_mask[$n_mask]=$wave[$i];
	$flux_mask[$n_mask]=$flux[$i];
	$n_mask++;
    }
}

$out_flux_pdl = interpol(pdl(@wave), pdl(@wave_mask), pdl(@flux_mask));
@out_flux=list($out_flux_pdl);


open(OUT,">$out_spec");
for ($i=0;$i<$n;$i++) {
    $k=$i+1;
    for ($j=0;$j<$wave_column;$j++) {
	print OUT "$k ";
#	print "$k ";
    }
    print OUT "$wave[$i] $out_flux[$i]\n";
#    print "$wave[$i] $out_flux[$i] $flux[$i]\n";
}
close(OUT);
exit;
