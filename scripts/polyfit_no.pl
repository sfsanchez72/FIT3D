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
    print "USE: poly.pl INPUT1.FITS NPOLY OUTPUT.FITS [NX1 NX2]\n";
    exit;
}

$infile=$ARGV[0];
$npoly=$ARGV[1];
$outfile=$ARGV[2];
$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
$nx1=0;
$nx2=0;
if ($#ARGV>2) {
    $nx1=$ARGV[3];
    $nx2=$ARGV[4];
} else {
    $nx2=$nx;
}
 
$a_out=zeroes($nx,$ny);
for ($j=0;$j<$ny;$j++) {
    my @a_spec;
    my $w_out;
    $nt=0;
    $min=1e12;
    $max=-1e12;

    for ($i=0;$i<$nx;$i++) {
	$w_out[$i]=$i;
	if (($i>$nx1)&&($i<$nx2)) {
	    $w_in[$nt]=$i;
	    $a_spec[$nt]=$a_in->at($i,$j);
	    if ($min>$a_spec[$nt]) {
		$min=$a_spec[$nt];
	    }
	    if ($max<$a_spec[$nt]) {
		$max=$a_spec[$nt];
	    }
	    
	    $nt++;
	}
    }
    ($s_x,$coeff) = fitpoly1d(pdl(@w_in),pdl(@a_spec),$npoly);
    @c=list($coeff);
    for ($i=0;$i<$nx;$i++) {
	$out_spec[$i]=0;
	for ($k=0;$k<$npoly;$k++) {
	    $out_spec[$i]=$out_spec[$i]+$c[$k]*($i**$k);
#	    print "$k $c[$k] $out_spec[$i]\n";
	}
#	if ($min>$out_spec[$i]) {
#       $min=$out_spec[$i];
#	}
#	if ($max<$out_spec[$i]) {
#	    $max=$out_spec[$i];
#	}

	set($a_out,$i,$j,$out_spec[$i]);
    }

    
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
