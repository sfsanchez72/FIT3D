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
    print "USE: create_mask.pl FIBERFLAT.fits mask.txt percentage(0-1) [[NX_min] [NX_max]]\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$f=$ARGV[2];

#$fiber_flat=$ARGV[2];


$array=rfits($infile);
($nx,$ny)=$array->dims;
if ($#ARGV==4) {
    $nx_min=$ARGV[3];
    $nx_max=$ARGV[4];
} else {
    $nx_min=0;
    $nx_max=$nx;
}


for ($j=0;$j<$ny;$j++) {
    my $spec=$array->slice("$nx_min:$nx_max,$j");
    ($mean,$sig,$median,$min,$max,$adev,$rms) = stats($spec);
    $med[$j]=$median;
#    print "$j $med[$j]\n";
}
$pdl=pdl(@med);
($mean,$sig,$median,$min,$max,$adev,$rms) = stats($pdl);
#print "($mean,$sig,$median,$min,$max,$adev,$rms)\n";
$n_bad=0;
open(OUT,">$out_file");
for ($j=0;$j<$ny;$j++) {
    if (abs(($med[$j]-$median)/$median)<$f) {
	$out_array=1;
    } else {
	$out_array=0;
	$n_bad++;
    }
#    if (abs(($med[$j]-$median))<($f*$sig)) {
#	$out_array=1;
#    } else {
#	$out_array=0;
#	$n_bad++;
#    }
	
    $k=$j+1;
    print OUT "$k $out_array\n";
}
close(OUT);
print "$out_file written. $n_bad bad fibers of $ny\n";


exit;
