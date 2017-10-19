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

#@ThAr=(27,55,83,111,139,165,166,194,222,250,278,306);

if ($#ARGV<3) {
    print "USE: clean_cr.pl INPUT.ms.fits OUTPUT.ms.fits NBOX NSIGMA\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$nbox=$ARGV[2];
$nsigma=$ARGV[3];
$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
$a_out=$a_in;
$a_out->sethdr($h);

$nc_tot=0;
for ($i=0;$i<$nx;$i++) {
    $nc_now=0;
    for ($jj=0;$jj<$ny;$jj++) {
#	$jj=$ThAr[$n]-1;
	$jmin=$jj-$nbox;
	$jmax=$jj+$nbox;
	if ($jmin<0) {
	    $jmin=0;
	    $jmax=$jmin+2*$nbox;
	}
	if ($jmax>=$ny) {
	    $jmax=$ny-1;
	    $jmin=$jmax-2*$nbox;
	}
#	print "[$jmin,$jmax] ";
	my @sec;
	$k=0;
	for ($j=$jmin;$j<$jmax;$j++) {
	    $sec[$k]=$a_in->at($i,$j);
	    $k++;
	}
	$m=median(@sec);
	$s=sigma(@sec);
#	print "$m $s\n";
	$val=$a_in->at($i,$jj);
	if (abs($val-$m)>$nsigma*$s) {
	    set($a_out,$i,$jj,$m);
	    $nc_tot++;
	    $nc_now++;
	}
    }
    print "$i/$nx $nc_now $nc_tot\n";
}

$a_out->wfits($outfile);
exit;

