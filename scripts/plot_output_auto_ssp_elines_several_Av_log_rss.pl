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

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: plot_output_auto_ssp_elines_several_Av_log_rss.pl INPUT_FILE.FITS NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
    exit;
}

$input=$ARGV[0];
$NY=$ARGV[1];
$dev=$ARGV[2];
$def=0;
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}





$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny,$nz)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL3"};
$cdelt=$pdl->hdr->{"CDELT3"};
$crpix=$pdl->hdr->{"CRPIX3"};



#print "$crval $cdelt $crpix\n";
if ($cdelt==0) {
    $cdelt=1;
}
for ($k=0;$k<$nz;$k++) {
    $key="NAME".$k;
    $name[$k]=$pdl->hdr->{$key};
    for ($i=0;$i<$nx;$i++) {
	$wave[$i][$k]=$crval+$cdelt*($i+1-$crpix);
	$flux[$i][$k]=$pdl->at($i,$NY,$k);
	if ($flux[$i][$k]>$y_max) {
	    $y_max=$flux[$i][$k];
	}
	if ($flux[$i][$k]<$y_min) {
	    $y_min=$flux[$i][$k];
	}
#    print "$crval $cdelt $crpix $wave[$i] $flux[$i]\n";
    }
}


$wmin=$wave[0][0];
$wmax=$wave[$nx-1][0];
print "$wmin $wmax\n";
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}

if ($#ARGV==6) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $def=1;
}
$ref_def=0;
if ($#ARGV==8) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $ref_lines=$ARGV[7];
    $cut_lines=$ARGV[8];
    $ref_def=1;
    $def=1;
}

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

#print "$y_min $y_max\n";

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
$color=1;
for ($k=0;$k<$nz;$k++) {
    for ($i=0;$i<$nx;$i++) {
	$wave_now[$i]=$wave[$i][$k];
	$flux_now[$i]=$flux[$i][$k];
#	if ($k==4) {
#	    print "$wave_now[$i] $flux_now[$i]\n";
#	}
    }
    pgsci($color);
    if ($color==1) {
	pgslw(5);
    } else {
	pgslw(2);
    }
    pgline($nx,\@wave_now,\@flux_now);
    pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.05*($k+1)*($y_max-$y_min),0,0,$name[$k]);
    $color++;
    if ($color>15) {
	$color=1;
    }
#    pgsch(0.9);
#    pgsci(2);
}

pgclose;
pgend;

my @sec;
$nn=0;
for ($i=0;$i<$nx;$i++) {
    if (($wave[$i]>$wmin)&&($wave[$i]<$wmax)) {
	$sec[$nn]=$flux[$i];
	$nn++;
    }
}


$median=median(@sec);
$mean=mean(@sec);
$sig=sigma(@sec);

$flux=$mean*$n*($wave[1]-$wave[0]);

print "$median +- $sig (F=$flux)\n";
exit;
