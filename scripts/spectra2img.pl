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


if ($#ARGV<1) {
    print "USE: spec2img.pl list_of_spectra.txt image.fits\n";
    exit;
}

$in_file=$ARGV[0];
$out_file=$ARGV[1];
#$crval=$ARGV[3];
#$cdelt=$ARGV[4];

$nf=0;
open(FH,"<$in_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$file[$nf]=$line;
	@data=split(/\//,$line);
	$name[$nf]=$data[$#data];
	$nf++;
    }
}
close(FH);

for ($j=0;$j<$nf;$j++) {
    open(FH,"<$file[$j]");
    $nx=0;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($j==0) {
	    if ($nx==0) {
		$crval=$data[1];
	    }
	    if ($nx==1) {
		$cdelt=$data[1]-$crval;
	    }
	}
	$id[$nx]=$data[0];
	$w[$nx]=$data[1];
	$f[$nx]=$data[2];
	$nx++;
    }
    close(FH);
    if ($j==0) {
	$pdl_out=zeroes($nx,$nf);
    }
    for ($i=0;$i<$nx;$i++) {
	set($pdl_out,$i,$j,$f[$i]);
    }
    $head="NAME".$j;
    $$h{$head}=$name[$j];
}
$$h{"CRPIX1"}=1;
$$h{"CRVAL1"}=$crval;
$$h{"CDELT1"}=$cdelt;






$pdl_out->sethdr( $h );
$pdl_out->wfits($out_file);

print "DONE\n";

exit;

