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
    print "USE: spec2D_plot.pl INPUT_FILE.FITS nplot DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
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
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};
if ($cdelt==0) {
    $cdelt=1;
}
for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($i,0);
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
#    print "$wave[$i] $flux[$i]\n";
}


$wmin=$wave[0];
$wmax=$wave[$n-1];
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
    if ($min!=$max) {
	$y_min=$min;    
	$y_max=$max;
    }
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

pgbegin(0,$dev,0,0);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(5/$NY);           # Set character height
pgsubp(1,$NY);
for ($i=0;$i<$NY;$i++) {
    $W1=$wmin+$i*($wmax-$wmin)/$NY;
    $W2=$wmin+($i+1)*($wmax-$wmin)/$NY;

    $med=median(@flux);
    $sig=sigma(@flux);
#    $y_min=$med-0.5*$sig;
#    $y_max=$med+1.5*$sig;


    pgenv($W1,$W2,$y_min,$y_max,0,0);   
#    pglabel("Wavelength","Flux","");
    pgsci(2);
    pgline($nx,\@wave,\@flux);
#    pgsch(1.2);
    pgsch(5);
    pgsci(8);
    pgptxt($W1+0.02*($W2-$W1),$y_max-0.2*($y_max-$y_min),0,0,"$W1-$W2");
    pgsch(1.2);
    pgsci(1);
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
