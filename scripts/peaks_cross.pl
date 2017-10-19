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


if ($#ARGV<7) {
    print "USE: peaks_cross.pl RAW.fits Spectral_axis[0/1] PEAKS_FILE x_width y_width plot nplot y_shift_limit [y_shift_ini] [stress] [nsec_start] [nsec_width]\n";
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
#$nsearch=$ARGV[7];
#$trace_file=$ARGV[8];
$y_shift_limit=abs($ARGV[7]);
if ($#ARGV==8) {
    $y_shift=$ARGV[8];
}

$stress=1;
if ($#ARGV==9) {
    $y_shift=$ARGV[8];
    $stress=$ARGV[9];
}

$nsec_start=0;
$nsec_width=0;
if ($#ARGV==11) {
    $y_shift=$ARGV[8];
    $stress=$ARGV[9];
    $nsec_start=$ARGV[10];
    $nsec_width=$ARGV[11];
}
#$dmin=$ARGV[7];
#$imin=$ARGV[8];
$peaks_file_out="out_".$peaks_file;

print "Reading $peaks_file...";
$npeaks=0;
open(FH,"<$peaks_file");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$id[$npeaks]=$data[0];
	$peak_y_pixel[$npeaks]=$data[1];
	$peak_y_pixel_ini[$npeaks]=$data[1]+$y_shift;
	$peak_y_max[$npeaks]=$data[2]+$y_shift;
    	$npeaks++;	
#    } else {
#	print OUT "$line";
    }
}
close(FH);

for ($j=0;$j<$npeaks;$j++) {
    $old_peak_y_max[$j]=$peak_y_max[$j];
}
print "\n";
open(OUT,">$peaks_file_out");
for ($j=0;$j<$npeaks;$j++) {
    if ($j==0) {
	print "$j $old_peak_y_max[$j] $peak_y_max[$j]\n";
    } else {
	if ($stress!=1) {
	    $d=$old_peak_y_max[$j]-$old_peak_y_max[$j-1];
	    $peak_y_max[$j]=$peak_y_max[$j-1]+($d)*$stress;#+$y_shift;
#	    print "$j $old_peak_y_max[$j] $peak_y_max[$j]=$peak_y_max[$j-1]+($d)*$stress \n";
	}
    }
    print OUT "$id[$j] $peak_y_pixel[$j] $peak_y_pixel_ini[$j] $peak_y_max[$j]\n";
}
print "STRESS=$stress\n";
close(OUT);
print "DONE\n";
print "NPEAKS=$npeaks\n";


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
@tmp_array=read_img($infile);
print "Done\n";
if ($spec_axis==0) {
    @in_array=@tmp_array;
    $nx=$naxis[0];
    $ny=$naxis[1];
} else {
    for ($j=0;$j<$naxis[1];$j++) {
	for ($i=0;$i<$naxis[0];$i++) {
	    $in_array[$i][$j]=$tmp_array[$j][$i];
	}
    }
    $nx=$naxis[1];
    $ny=$naxis[0];
}

#
# We select the central spectral regions
#
$mean_point=int($nx/2);
#
# We cross-correlate to get the $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #d shift
#
print "Starting Cross-Correlation\n";
$y_shift_tot=$y_shift;
$nloop=0;
do {
$mm=$mean_point;
my @sec_array;
    $min=1e12;
    $max=-1e12;

for ($j=0;$j<$ny;$j++) {
    $sec_array[$j]=0;
    $peaks_array[$j]=0;
    if ($width>0) {
	for ($i=($mm-$width);$i<($mm+$width);$i++) {
	    $sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	}
    } else {
	for ($i=$mm;$i<($mm+1);$i++) {
	    $sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
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
for ($k=0;$k<$npeaks;$k++) {
    $centroid=$peak_y_max[$k];
    $peak_flux[$k]=$sec_array[int($centroid)];
#    print "$k $peak_flux[$k]\n";
}
$median=median(@peak_flux);
if ($max>2*$median) {
    $max=2*$median;
}
#print "$median\n";

$sigma=($y_width/(2*2.345));
for ($k=0;$k<$npeaks;$k++) {
    $centroid=$peak_y_max[$k];
    for ($j=0;$j<$ny;$j++) {
#	$peaks_array[$j]=$peaks_array[$j]+$sec_array[int($centroid)]*exp(-0.5*(($j-$centroid)/$sigma)**2);	
	$peaks_array[$j]=$peaks_array[$j]+$median*exp(-0.5*(($j-$centroid)/$sigma)**2);	
    }
#    print "$cut[$j] $sec_array[$j] $min $max $peaks_array[$j]\n";
}

if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsubp(1,$nplot);
	for ($i=0;$i<$nplot;$i++) {
	    pgenv(($ny/$nplot)*$i,($ny/$nplot)*($i+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","");
	    pgline($ny,\@cut,\@sec_array);    
	    pgsci(3);
	    pgline($ny,\@cut,\@peaks_array);    
	    pgsch(4.0);           # Set character height
	    for ($k=0;$k<$npeaks;$k++) {
		pgsci(7);
		$x=$peak_y_pixel[$k];
		$y=0.8*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],1);
		pgsci(2);
		$x=$peak_y_pixel_new[$k][$mm];
		$y=0.7*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],5);
		pgsci(4);
		$x=$peak_y_max[$k];
		$y=0.6*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],2);
	    }
	    pgsci(3);
	    pgline(2,[0,$ny],[$imin*$max,$imin*$max]);
	    pgsci(1);
	    pgsch(1.2);           # Set character height
	}
	pgsci(1);
	pgclose;
	pgend;
	if ($command ne "A") {
	    print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	    chop($command);
	}
	if ($command eq "0") {
	    $plot=0;
	}
    }

$k=20;
$n_max=2**$k;
while ($n_max>$ny) {
    $k--;
    $n_max=2**$k;
#    print "$k $n_max $ny\n";
}
#$n_max=512;
if ($nsec_width==0) {
    $nsec_width=$n_max;
}
for ($j=0;$j<$nsec_width;$j++) {
    $c_sec_array[$j]=$sec_array[$j+$nsec_start];
    $c_peaks_array[$j]=$peaks_array[$j+$nsec_start];
}

my $fft1=new Math::FFT(\@c_sec_array);
my $fft2=new Math::FFT(\@c_peaks_array);
my $corr= $fft1->correl($fft2);
@a_tmp=@$corr;
$max=-10000;
$min=1e17;
$sum=0;
for ($j=0;$j<$nsec_width;$j++) {    
    if ($j<($nsec_width/2)) {
	$a_corr[$j]=abs($a_tmp[$j+int($nsec_width/2)]);
    } else {
	$a_corr[$j]=abs($a_tmp[$j-int($nsec_width/2)]);
    }
#    $a_corr[$j]=$a_tmp[$j];
	$a_xx[$j]=$j;
#	print "$a_corr[$j] $a_xx[$j]\n";
	if ($a_corr[$j]<$min) {
	    $min=$a_corr[$j];
	}
	if ($a_corr[$j]>$max) {
	    $max=$a_corr[$j];
	}
}
$max=-1e12;
for ($j=($nsec_width/2-$y_width*2);$j<($nsec_width/2+$y_width*2);$j++) {
    if ($max<$a_corr[$j]) {
	$max=$a_corr[$j];
	$peak_point=$j;
    }   	
    $a=$peak_point-1;
    $b=$peak_point;
    $c=$peak_point+1;
    $fa=-$a_corr[$a];
    $fb=-$a_corr[$b];
    $fc=-$a_corr[$c];
    $den=($fc-2*$fb+$fa);
    if ($den!=0) {
	$real_peak_point=$c-($b-$a)*(($fc-$fb)/$den+0.5);
    } else {
	$real_peak_point=$j;
    }
}
$y_shift=$real_peak_point-int($nsec_width/2);
print "Y_shift=$y_shift ($y_shift_limit)\n";

if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($nsec_width/2-$y_width*2,$nsec_width/2+$y_width*2,$min,$max,0,0);
    pgline($nsec_width,\@a_xx,\@a_corr);    
    pgsci(1);
    pgclose;
    pgend;
    if ($command ne "A") {
	print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	chop($command);
    }
    if ($command eq "0") {
	$plot=0;
    }
}

for ($k=0;$k<$npeaks;$k++) {
    $peak_y_pixel_ini[$k]=$peak_y_pixel_ini[$k]+$y_shift;
    $peak_y_max[$k]=$peak_y_max[$k]+$y_shift;
}
$y_shift_tot=$y_shift_tot+$y_shift;
$nloop++;
} while ((abs($y_shift)>$y_shift_limit)&&($nloop<10));
print "Total Shift=$y_shift_tot\n";
print "DONE\n";

open(FH,">>peak_cross.out"); 
 print FH "$infile $peaks_file $y_shift_tot\n";
close(FH);

exit;



