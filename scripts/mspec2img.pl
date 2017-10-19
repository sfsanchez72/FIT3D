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


if ($#ARGV<1) {
    print "USE: spec2img.pl spec_list.txt image.fits [crval] [cdelt]\n";
    exit;
}

$spec_list=$ARGV[0];
$out_file=$ARGV[1];
$def=0;
if ($#ARGV==3) {
    $crval=$ARGV[2];
    $cdelt=$ARGV[3];
    $def=1;
}

$ny=0;
open(LIST,"<$spec_list");
while($line=<LIST>) {
    chop($line);
    $file[$ny]=$line;
    $ny++;
}
close(LIST);

for ($j=0;$j<$ny;$j++) {
    $n_unc=0;
    $unc_file=$file[$j];
    open(FH,"<$unc_file");
    while($line=<FH>) {
	@data=split(" ",$line);
	$index_unc[$n_unc]=$data[0];
	$wave_unc[$n_unc]=$data[1];
	$out_array[$j][$n_unc]=$data[2];
#	$flux_unc_ini[$n_unc]=$data[2];
#    print "$wave_unc[$n_unc] $flux_unc[$n_unc]\n";
	$n_unc++;
    }
    close(FH);
}
if ($def==0) {
    $crval=$wave_unc[0];
    $cdelt=$wave_unc[1]-$wave_unc[0];
}


print "Writting the 2D image $out_file\n";
system("rm $out_file");
write_rss($out_file,$n_unc,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;

