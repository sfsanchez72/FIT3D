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


if ($#ARGV<8) {
    print "USE: peak_find_width.pl RAW.fits Spectral_axis[0/1] Coadd_width plot nplot nsearch DMIN IMIN(% of the MAX) OUTFILE.WIDTH [colum_search]\n";
    exit;
}


$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$width=$ARGV[2];
$plot=$ARGV[3];
$nplot=$ARGV[4];
$nsearch=$ARGV[5];
$dmin=$ARGV[6];
$imin=$ARGV[7];
#print "$imin\n";
$outfile=$ARGV[8];

#if ($nsearch==2*int($nsearch/2)) {
#    $nsearch++;
#}

#print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
@tmp_array=read_img($infile);
#print "Done\n";

#$a_in=rfits($infile);
#$h=$a_in->gethdr;
#@naxis=$a_in->dims;
#@tmp_array=list($a_in);

if ($spec_axis==0) {
#    @in_array=list($a_in);
    @in_array=@tmp_array;
    $nx=$naxis[0];
    $ny=$naxis[1];
} else {
    for ($j=0;$j<$naxis[1];$j++) {
	for ($i=0;$i<$naxis[0];$i++) {
	    $in_array[$i][$j]=$tmp_array[$j][$i];
#	    $in_array[$i][$j]=$a_in->at($i,$j);#tmp_array[$j][$i];
	}
    }
    $nx=$naxis[1];
    $ny=$naxis[0];
}

#
# We select the central spectral regions
#
$mean_point=int($nx/2);
$min=1e12;
$max=-1e12;
for ($j=0;$j<$ny;$j++) {
    $sec_array[$j]=0;
    for ($i=($mean_point-$width);$i<($mean_point+$width);$i++) {
	$sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
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
# We look for the peaks
#
#print "Looking for peaks with in $nsearch pixels...";
$npeaks=0;
for ($j=$nsearch;$j<($ny-$nsearch);$j++) {
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
	$rrr=$imin*$max;
	if ($sec_array[$j]<($imin*$max)) {
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

	    
#	$a=$j-1;
#	$b=$j;
#	$c=$j+1;
#	$fa=$sec_array[$a];
#	$fb=$sec_array[$b];
#	$fc=$sec_array[$c];
#	$peak_y_max[$npeaks]=$b-0.5*((($b-$a)**2)*($fa-$fc)-(($b-$c)**2)*($fb-$fa))/(($b-$a)*($fb-$fc)-($b-$c)*($fb-$fa));

	$a=$j-1;
	$b=$j;
	$c=$j+1;
	$fa=-$sec_array[$a];
	$fb=-$sec_array[$b];
	$fc=-$sec_array[$c];
	$den=($fc-2*$fb+$fa);
	if ($den!=0) {
	    $peak_y_max[$npeaks]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	} else {
	    $peak_y_max[$npeaks]=0;
	}
#	print "$peak_y_max[$npeaks]=$c-($b-$a)*(($fc-$fb)/$den+0.5)\n";
	$npeaks++;
    }
	    
}
#print "DONE\n";
#print "NPEAKS=$npeaks\n";
#
# We save the result
#



open(CON,">peak_find_width.config");
print CON "0 $npeaks 0.01 0.001\n";
for ($i=0;$i<$npeaks;$i++) {
    $j=$peak_y_max[$i];
    print CON "eline\n";
    print CON "$j       0       0       0       -1\n";
    print CON "1e6     1       -0.1    1e30    -1\n";
#    if ($i==0) {
	print CON "1        1       0.1     4      -1\n";
#    } else {
#	print CON "1       1       0        0       1\n";
#    }
    print CON "0        0       0       0       -1\n";
    print CON "0        0       0       0       -1\n";
    print CON "0        0       0       0       -1\n";
    print CON "0        0       0       0       -1\n";
    print CON "0        0       0       0       -1\n";
    print CON "0        0       0       0       -1\n";
}
close(CON);

$nbox=3;
open(FH,">sec_array.spec");
for ($i=0;$i<$ny;$i++) {
    print FH "$i $cut[$i] $sec_array[$i]\n";
    $mask[$i]=0;
    for ($j=0;$j<$npeaks;$j++) {
	$d=abs($peak_y_max[$j]-$cut[$i]);
	if ($d<$nbox) {
	    $mask[$i]=1;
	}
    }

}
close(FH);

$m1=0;
$m2=0;

open(FH,">sec_array.mask");
for ($i=1;$i<$ny;$i++) {
    if (($mask[$i]==1)&&($mask[$i-1]==0)) {
	$m2=$cut[$i-1];
	print FH "$m1 $m2\n";
    }
    if (($mask[$i]==0)&&($mask[$i-1]==1)) {
	$m1=$cut[$i];
    }

}
close(FH);

#
# We create the mask
#



if ($plot==1) {
    $dev="/xs";
} else {
    if ($plot==0) {
	$dev="/null";
    } else {
	$dev="peak_find_width.ps/CPS";
    }
}




$call="fit_spec_back.pl sec_array.spec none 0 10000 none 0 sec_array.mask peak_find_width.config ".$dev." > fit.junk";
system($call);

$k=0;
open(OUT,"<out.fit_spectra");
open(FH,">$outfile");
print FH "# PEAKS in file $in_file\n";
print FH "# ID pix_peak peak sigma\n";
while($out=<OUT>) {
    chop($out);
    if ($out =~ "eline") {
	@data=split(" ",$out);
	$sigma[$k]=$data[5];
	print FH "$k $peak_y_pixel[$k] $peak_y_max[$k] $sigma[$k]\n";
	$k++;
    }
}
close(OUT);
$mean_sigma=mean(@sigma);
$sigma_sigma=sigma(@sigma);
$fwhm=2.345*$mean_sigma;
$s_fwhm=2.345*$sigma_sigma;
print "$npeaks $mean_sigma $sigma_sigma $fwhm $s_fwhm\n";
#print FH "-1 -1 -1 $mean_sigma $sigma_sigma $fwhm\n";
close(FH);

$N=0;
open(FH,"<out_mod_res.fit_spectra");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $X[$N]=$data[0];
    $Y[$N]=$data[2];
    $N++;
}
close(FH);


#
# We plot the section
#

    pgbegin(0,$dev,1,1);
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
	    pgsci(2);
	    $x=$peak_y_pixel[$k];
	    $y=0.8*$max;
	    pgpoint(1,[$x],[$y],5);
	    pgsci(4);
	    $x=$peak_y_max[$k];
	    $y=0.6*$max;
	    pgpoint(1,[$x],[$y],2);
	}
	pgsci(3);
	pgline(2,[0,$ny],[$imin*$max,$imin*$max]);
	pgsci(8);
	pgline($N,\@X,\@Y);    
	pgsci(1);
	pgsch(1.2);           # Set character height
    }
    pgsci(1);
    pgclose;
    pgend;
if ($plot==1) {
    print "Press Enter"; <stdin>;
}


exit;



