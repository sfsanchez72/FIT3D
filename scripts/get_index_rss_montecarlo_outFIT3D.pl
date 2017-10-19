#!/usr/bin/perl

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

#use POSIX;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<4) {
    print "USE:get_index_rss_montecarlo.pl SPEC.RSS.fits NOISE.RSS.fits NSIM auto_ssp.rss.out DEV\n";    
    exit;
}

$spec_file=$ARGV[0];
$noise_file=$ARGV[1];
$NSIM=$ARGV[2];
$file_z=$ARGV[3];
$dev=$ARGV[4];

$n=0;
open(FH,"<$file_z");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$z[$n]=$data[4];
	$n++;	
    }
}
close(FH);


$a=rfits($spec_file);
($nx,$ny)=$a->dims;
for ($j=0;$j<$ny;$j++) {
    $call="img2spec.pl ".$spec_file." ".$j." get_index_rss.spec";
    system($call);
    $call="img2spec.pl ".$noise_file." ".$j." get_index_noise.spec";
    system($call);
    $call="get_index_montecarlo.pl get_index_rss.spec get_index_noise.spec ".$NSIM." ".$z[$j]." ".$dev;
    system($call);
}


exit;
