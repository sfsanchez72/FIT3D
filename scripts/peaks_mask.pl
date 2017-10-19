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


if ($#ARGV<9) {
    print "USE: peaks_mask.pl RAW.fits Spectral_axis[0/1] PEAKS_FILE x_width y_width plot nplot nsearch mask.txt %FLUX_limit [y_shift] [plot_width]\n";
    exit;
}

$y_shift=0;
$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$peaks_file=$ARGV[2];
$width=$ARGV[3];
if ($width==0) {
    $width=1;
}
$y_width=$ARGV[4];
$plot=$ARGV[5];
$nplot=$ARGV[6];
$nsearch=$ARGV[7];
$trace_file=$ARGV[8];
$per_limit=$ARGV[9];
if ($#ARGV==10) {
    $y_shift=$ARGV[10];
}
if ($#ARGV==11) {
    $y_shift=$ARGV[10];
    $plot_width=$ARGV[11];
}
#$dmin=$ARGV[7];
#$imin=$ARGV[8];

print "Reading $peaks_file...";
$npeaks=0;
open(FH,"<$peaks_file");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$id[$npeaks]=$data[0];
	$peak_y_pixel[$npeaks]=$data[1];
#	$peak_y_pixel_ini[$npeaks]=$data[1]+$y_shift;
	$peak_y_max[$npeaks]=$data[2]+$y_shift;
	$peak_y_pixel_ini[$npeaks]=int($peak_y_max[$npeaks]);
	$npeaks++;
    }
}
close(FH);
print "DONE\n";
print "NPEAKS=$npeaks\n";


print "Reading file $infile ";
#$nax=read_naxes($infile);   

$tmp=rfits($infile);
@naxis=$tmp->dims();

#@tmp_array=read_img($infile);
#@tmp_array=list($pdl);
print "Done\n";

if ($spec_axis==0) {
#    @in_array=list($pdl);
    $nx=$naxis[0];
    $ny=$naxis[1];
    $pdl=$tmp;
} else {
#    for ($j=0;$j<$naxis[1];$j++) {
#	for ($i=0;$i<$naxis[0];$i++) {
#	    $in_array[$i][$j]=$pdl->at($j,$i);
#	}
#    }
#        $pdl=$tmp;
    $pdl=$tmp->xchg(0,1);
    $nx=$naxis[1];
    $ny=$naxis[0];
}

#print "spec_axis=$spec_axis\n";
#print "NN= @naxis $nx,$ny\n";

#
# We select the central spectral regions
#
$mean_point=int($nx/2);
print "We start to mask\n";
#
# 1st Half
#
@peak_y_pixel=@peak_y_pixel_ini;
$mm=$mean_point;
#for ($mm=$mean_point;$mm<($nx-$width);$mm++) {
#    print "COLUMN=$mm\n";
$min=1e12;
$max=-1e12;
$nmin=$mm-$width;
$nmax=$mm+$width;
$nnn=$nmax-$nmin+1;
my @sec_array;
$sec_pdl=$pdl->slice("$nmin:$nmax,:");
@dims=$sec_pdl->dims;
@dims_old=$pdl->dims;
#print "DIM=@dims @dims_old\n";

#exit;
for ($j=0;$j<$ny;$j++) {
    $sec_array[$j]=0;
	if ($width>0) {
	    for ($i=0;$i<$nnn;$i++) {
#		print "$i,$j\n";
		$sec_array[$j]=$sec_array[$j]+$sec_pdl->at($i,$j);
	    }
	} else {
	    for ($i=0;$i<1;$i++) {
		$sec_array[$j]=$sec_array[$j]+$sec_pdl->at($i,$j);
	    }
	}
	if ($min>$sec_array[$j]) {
	    $min=$sec_array[$j];
	}
	if ($max<$sec_array[$j]) {
	    $max=$sec_array[$j];
	}	    
	$cut[$j]=$j;
    }
#
# We look for the peaks around the original ones...
#
#print "Looking for peaks with in $nsearch pixels...";

    for ($k=0;$k<$npeaks;$k++) {
	$existe_peak=0;
	    if ($mm==$mean_point) {
		$peak_y_pre=$peak_y_pixel[$k];
	    } else {
		$peak_y_pre=$peak_y_pixel_new[$k][$mm-1];
	    }
#	for ($j=($peak_y_pixel[$k]-$y_width);$j<($peak_y_pixel[$k]+$y_width);$j++) {
	for ($j=($peak_y_pre-$y_width);$j<($peak_y_pre+$y_width);$j++) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($sec_array[$j-$i]<$sec_array[$j-$i-1]) {
		    $peak=0;
		}
		if ($sec_array[$j+$i]<$sec_array[$j+$i+1]) {
		    $peak=0;
		}	       
	    }
	    if ($peak==1) {
		$peak_y_pixel_new[$k][$mm]=$j*1.0;		
		$existe_peak=1;
	    } 

	}
	if ($existe_peak==0) {
	    if ($mm==$mean_point) {
		$peak_y_pixel_new[$k][$mm]=$peak_y_pixel[$k]*1.0;				$peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	    } else {
		$peak_y_pixel_new[$k][$mm]=$peak_y_pixel_new[$k][$mm-1];
		$peak_y_max_new[$k][$mm]=$peak_y_max_new[$k][$mm-1];
	    }
	} else {
	    $j=$peak_y_pixel_new[$k][$mm];		
	    $a=$j-1;
	    $b=$j;
	    $c=$j+1;
	    $fa=-$sec_array[$a];
	    $fb=-$sec_array[$b];
	    $fc=-$sec_array[$c];
	    $den=($fc-2*$fb+$fa);
	    if ($den!=0) {
#		$peak_y_max_new[$k][$mm]=$a-0.5*((($b-$a)**2)*($fa-$fc)-(($b-$c)**2)*($fb-$fa))/(($b-$a)*($fb-$fc)-($b-$c)*($fb-$fa));	    
		$peak_y_max_new[$k][$mm]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	    } else {
		$peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	    }
	}
	$peak_y_pixel[$k]=$peak_y_pixel_new[$k][$mm];
	$flux_peak[$k]=$fb*(-1);
#	print "$k $fb \n";
    }
$median_flux_peak=median(@flux_peak);
$low=$median_flux_peak*($per_limit);
$min=0;
$max=3*$median_flux_peak;
#print "$median_flux_peak $low\n";
#
# We plot the section
#
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgask(0);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgsubp(1,$nplot);
	for ($i=0;$i<$nplot;$i++) {
	    pgenv(($ny/$nplot)*$i,($ny/$nplot)*($i+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","");
	    pgline($ny,\@cut,\@sec_array);    
	    pgsch(4.0);           # Set character height
	    for ($k=0;$k<$npeaks;$k++) {
		pgsci(7);
#		$x=$peak_y_pixel[$k];
		$x=$peak_y_max[$k];
		$y=0.8*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],1);
		pgsci(2);
		$x=$peak_y_pixel_new[$k][$mm];
		$y=0.7*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],5);
		pgsci(4);
		$x=$peak_y_max_new[$k][$mm];
		$y=0.6*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],2);
	    }
	    pgsci(3);
	    pgline(2,[0,$ny],[$low,$low]);
	    pgsci(1);
	    pgsch(1.2);           # Set character height
	}
	pgsci(1);
	pgclose;
	pgend;
}

$nmasked=0;

open(FH,">$trace_file");
for ($k=0;$k<$npeaks;$k++) {
    $mm=$k+1;
    if ($flux_peak[$k]>$low) {
	$mask=1;
    } else {
	$nmasked++;
	$mask=0;
    }
    print FH "$mm $mask\n";
}
close(FH);

print "$nmasked peaks\n";

exit;
