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
    print "USE: spec2img.pl spec.txt image.fits ny crval cdelt\n";
    exit;
}

$unc_file=$ARGV[0];
$out_file=$ARGV[1];
$ny=$ARGV[2];
$crval=$ARGV[3];
$cdelt=$ARGV[4];

$n_unc=0;
open(FH,"<$unc_file");
while($line=<FH>) {
    @data=split(" ",$line);
    $index_unc[$n_unc]=$data[0];
    $wave_unc[$n_unc]=$data[1];
    $flux_unc_ini[$n_unc]=$data[2];
#    print "$wave_unc[$n_unc] $flux_unc[$n_unc]\n";
    $n_unc++;
}
close(FH);

for ($i=0;$i<$n_unc;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$out_array[$j][$i]=$flux_unc_ini[$i];
    }
}

print "Writting the 2D image $out_file\n";
system("rm $out_file");
write_rss($out_file,$n_unc,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;

