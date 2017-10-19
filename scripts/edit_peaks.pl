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


if ($#ARGV<4) {
    print "USE: edit_peaks.pl RAW.fits Spectral_axis[0/1] INPUT_PEAKS_FILE OUTPUT_PEAKS_FILE y_width\n";
    exit;
}

$y_shift=0;
$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$peaks_file=$ARGV[2];
$out_peaks_file=$ARGV[3];
$y_width=$ARGV[4];

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
    }
}
close(FH);
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
pgbegin(0,"/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv(0,$ny,$min,$max,0,0);
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
pgsci(1);
pgsch(1.2);           # Set character height
pgsci(1);
pgclose;
pgend;
print "Press Enter"; <stdin>;


exit;
