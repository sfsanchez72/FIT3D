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


if ($#ARGV<4) {
    print "USE: spectra2img_norm.pl list_of_spectra.txt image.fits wave_norm crval cdelt [F_NORM]\n";
    exit;
}

$in_file=$ARGV[0];
$out_file=$ARGV[1];
$wn=$ARGV[2];
$crval=$ARGV[3];
$cdelt=$ARGV[4];
$factor=1;
if ($#ARGV==5) {
    $factor=$ARGV[5];
}

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
    $dist=1e12;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($#data>1) {
            $id[$nx]=$data[0];
            $w[$nx]=$data[1];
            $f[$nx]=$data[2];
        } else {
            $w[$nx]=$data[0];
            $f[$nx]=$data[1];
        }
	if (abs($w[$nx]-$wn)<$dist) {
	    $FN=$f[$nx];
	    $dist=abs($w[$nx]-$wn);
	}
	$nx++;

    }
    close(FH);
    if ($j==0) {
	$NX=int(($w[$nx-1]-$crval)/$cdelt);
	for ($i=0;$i<$NX;$i++) {
	    $wave_new[$i]=$crval+$cdelt*$i;
	}
	$pdl_out=zeroes($NX,$nf);
    }
    $pdl_new=interpol(pdl(@wave_new),pdl(@w),pdl(@f));
    for ($i=0;$i<$NX;$i++) {
	$val=$pdl_new->at($i);
	$val=$val/$FN;
	set($pdl_out,$i,$j,$val);
    }
    $head="NAME".$j;
    $$h{$head}=$name[$j];
    $head="NORM".$j;
    $$h{$head}=$FN*$factor;
}
$$h{"WAVENORM"}=$wn;
$$h{"CRPIX1"}=1;
$$h{"CRVAL1"}=$crval;
$$h{"CDELT1"}=$cdelt;


$pdl_out->sethdr( $h );
$pdl_out->wfits($out_file);

print "DONE\n";

exit;

