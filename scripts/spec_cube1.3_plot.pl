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

if ($#ARGV<2) {
    print "USE: spec_cube_plot.pl INPUT_FILE.CUBE.FITS NX NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
    exit;
}

$input=$ARGV[0];
$NX=$ARGV[1];
$NY=$ARGV[2];
$dev=$ARGV[3];
$def=0;
if ($#ARGV==5) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $def=1;
}





$y_min=1e12;
$y_max=-1e12;
$n=0;
$input_now=$input."[0]";
$pdl=rfits("$input_now");
$input_now=$input."[1]";
$pdl_var=rfits("$input_now");
$input_now=$input."[3]";
$pdl_mask=rfits("$input_now");
($nz,$ny,$nx)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL3"};
$cdelt=$pdl->hdr->{"CDELT3"};
$crpix=$pdl->hdr->{"CRPIX3"};

#print "$crval $cdelt $crpix\n";
if ($cdelt==0) {
    $cdelt=1;
}
@stats=stats($pdl);

open(OUT,">spec_cube1.3_plot.out");
for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($NX,$NY,$i);
    $var[$i]=$pdl_var->at($NX,$NY,$i);
    $var1[$i]=$flux[$i]+$pdl_var->at($NX,$NY,$i);
    $var2[$i]=$flux[$i]-$pdl_var->at($NX,$NY,$i);
    $mask[$i]=($stats[0]->at(0))*$pdl_mask->at($NX,$NY,$i);
#    if ($mask[$i]>0) {
#	$flux[$i]=0;
#    } 
#    $mask[$i]=$pdl_mask->at($NX,$NY,$i);
    print OUT "$i $wave[$i] $flux[$i] $var[$i]\n";
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
#    print "$crval $cdelt $crpix $wave[$i] $flux[$i]\n";
}
close(OUT);

$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==5) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $def=1;
}

if ($#ARGV==7) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $wmin=$ARGV[6];
    $wmax=$ARGV[7];
    $def=1;
}
$ref_def=0;
if ($#ARGV==9) {
    $min=$ARGV[4];
    $max=$ARGV[5];
    $wmin=$ARGV[6];
    $wmax=$ARGV[7];
    $ref_lines=$ARGV[8];
    $cut_lines=$ARGV[9];
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
pgsls(4);
pgsci(8);
pgline($nx,\@wave,\@var1);
pgline($nx,\@wave,\@var2);
pgsci(2);
pgsls(1);
pgline($nx,\@wave,\@mask);
pgsci(1);
pgsls(1);
pgslw(3);
pgline($nx,\@wave,\@flux);
pgslw(1);
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
