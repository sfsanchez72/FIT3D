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

if ($#ARGV<3) {
    print "USE: spec_cube_plot.pl INPUT_FILE.CUBE.FITS INPUT_FILE2.CUBE.fits NX NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
    exit;
}

$input=$ARGV[0];
$input2=$ARGV[1];
$NX=$ARGV[2];
$NY=$ARGV[3];
$dev=$ARGV[4];
$def=0;
if ($#ARGV==6) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $def=1;
}





$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
$pdl2=rfits("$input2");
($nz,$ny,$nx)=$pdl->dims;
($nz2,$ny2,$nx2)=$pdl2->dims;
$crval=$pdl->hdr->{"CRVAL3"};
$cdelt=$pdl->hdr->{"CDELT3"};
$crpix=$pdl->hdr->{"CRPIX3"};

$crval2=$pdl2->hdr->{"CRVAL3"};
$cdelt2=$pdl2->hdr->{"CDELT3"};
$crpix2=$pdl2->hdr->{"CRPIX3"};

#print "$crval $cdelt $crpix\n";
if ($cdelt==0) {
    $cdelt=1;
}
for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($NX,$NY,$i);
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
}


for ($i=0;$i<$nx2;$i++) {
    $wave2[$i]=$crval2+$cdelt2*($i+1-$crpix2);
    $flux2[$i]=$pdl2->at($NX,$NY,$i);
#    print "$wave2[$i] $flux2[$i] $wave[$i] $flux[$i]\n";
}


$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==6) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $def=1;
}

if ($#ARGV==8) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $wmin=$ARGV[7];
    $wmax=$ARGV[8];
    $def=1;
}
$ref_def=0;
if ($#ARGV==10) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $wmin=$ARGV[7];
    $wmax=$ARGV[8];
    $ref_lines=$ARGV[9];
    $cut_lines=$ARGV[10];
    $ref_def=1;
    $def=1;
}

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

#print "$y_min $y_max\n";

if ($ref_def==1) {
    $nl=0;
    open(FH,"<$ref_lines");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($data[0]>$cut_lines) {
	    $wave_line[$nl]=$data[1];
	    if ($nl>0) {
		if (abs($wave_line[$nl]-$wave_line[$nl-1])>5) {
		    $nl++;
		    }
	    } else {
		$nl++;
	    }
	}
    }
    close(FH);
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
pgline($nx,\@wave,\@flux);
pgsci(2);
pgline($nx2,\@wave2,\@flux2);
pgsci(1);
pgsch(0.9);
pgsci(2);
for ($j=0;$j<$nl;$j++) {
    if (($wave_line[$j]>$wmin)&&($wave_line[$j]<$wmax)) {
	if ($j==(2*int($j/2))) {
	    pgptxt($wave_line[$j],0.8*$y_max,90,0,"$wave_line[$j]");
	} else {
	    pgptxt($wave_line[$j],0.7*$y_max,90,0,"$wave_line[$j]");
	}
    }
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
