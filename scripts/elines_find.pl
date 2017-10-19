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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<3) {
    print "USE: elines_find.pl spec.txt nsearch LIMIT(% max) OUTFILE [WMIN WMAX]\n";
    exit;
}

$input=$ARGV[0];
$nsearch=$ARGV[1];
$imin=$ARGV[2];
$outfile=$ARGV[3];
if ($#ARGV==5) {
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
} else {
    $wmin=0;
    $wmax=1e12;
}
$dmin=0;
$y_min=1e12;
$y_max=-1e12;
$n=0;
open(FH,"<$input");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $flux[$n]=$data[2];
    $wave[$n]=$data[1];
    if (($wave[$n]>$wmin)&&($wave[$n]<$wmax)) {

    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n];
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n];
    }
	$n++;
    }
}
close(FH);
$y_min=0;
$wmin=$wave[0];
$wmax=$wave[$n-1];
$crval=$wmin;
$cdelt=$wave[1]-$wave[0];


$npeaks=0;
for ($j=$nsearch;$j<($n-$nsearch);$j++) {
    $peak=1;
    for ($i=0;$i<$nsearch;$i++) {
	if ($flux[$j-$i]<$flux[$j-$i-1]) {
	    $peak=0;
	}
	if ($flux[$j+$i]<$flux[$j+$i+1]) {
	    $peak=0;
	}
    }
    if ($peak==1) {
	$rrr=$imin*$y_max;
	if ($flux[$j]<($imin*$y_max)) {
	    $peak=0;
	}
    }

    if ($peak==1) {
	if ($npeaks>0) {
	    $delta=$j-$peak_y_pixel[$npeaks-1];
	    if ($delta<$dmin) {
		$peak=0;
	    }
	}
    }

    if ($peak==1) {
	$peak_y_pixel[$npeaks]=$j;

	    
	$a=$j-1;
	$b=$j;
	$c=$j+1;
	$fa=-$flux[$a];
	$fb=-$flux[$b];
	$fc=-$flux[$c];
	$den=($fc-2*$fb+$fa);
	if ($den!=0) {
	    $peak_y_max[$npeaks]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	} else {
	    $peak_y_max[$npeaks]=0;
	}
	$npeaks++;
    }
	    
}
#print "DONE\n";
#print "NPEAKS=$npeaks\n";
#
# We save the result
#
open(FH,">$outfile");
print FH "# Elines in file $in_file\n";
for ($j=0;$j<$npeaks;$j++) {
    $wave_peak[$j]=$crval+$cdelt*$peak_y_max[$j];
    print FH "$wave_peak[$j]\n";
}
close(FH);


#
# We plot the section
#
#if ($plot==1) {
#    pgbegin(0,"/xs",1,1);
    pgbegin(0,"elines_find.ps/CPS",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
#    pgsubp(1,$nplot);
#    for ($i=0;$i<$nplot;$i++) {
    pgenv($wmin,$wmax,$y_min,$y_max,0,0);
    pglabel("Y-axis","Counts","");
    pgline($n,\@wave,\@flux);    
    pgsch(1.0);           # Set character height
    for ($k=0;$k<$npeaks;$k++) {
	pgsci(2);
	$x=$crval+$cdelt*$peak_y_pixel[$k];
	$y=0.8*$y_max;
	pgpoint(1,[$x],[$y],5);
	pgsci(4);
	$x=$wave_peak[$k];
	$y=0.6*$y_max;
	pgpoint(1,[$x],[$y],2);
    }
    pgsci(3);
    pgline(2,[$wmin,$wmax],[$imin*$y_max,$imin*$y_max]);
    pgsci(1);
    pgsch(1.2);           # Set character height
#    }
    pgsci(1);
    pgclose;
    pgend;
#    print "Press Enter"; <stdin>;
#}
exit;



