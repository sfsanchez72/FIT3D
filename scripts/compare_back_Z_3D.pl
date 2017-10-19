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
use PDL::Image2D;



#$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<10) {
    print "USE: compare_back_Z.pl RSS.fits PT.txt BACK_SPEC OUTFILE MASK_LIST VEL_MIN VEL_MAX DELTA_VEL SIGMA_DISP WAVE_SCALE PLOT [min max] [wmin wmax]\n";
    exit;
}

$unc_file=$ARGV[0];
$pt=$ARGV[1];
$c_file=$ARGV[2];
$outfile=$ARGV[3];
$mask_list=$ARGV[4];
$vel_min=$ARGV[5];
$vel_max=$ARGV[6];
$delta_vel=$ARGV[7];
$sigma=$ARGV[8];
$wave_scale=$ARGV[9];
$out_file="junk.junk";
$factor=1;
$box=1;
$plot=$ARGV[10];
$smooth=1;
$MIN_CHISQ=1e12;

$def=0;
if ($#ARGV==12) {
    $min=$ARGV[11];
    $max=$ARGV[12];
    $def=1;
}

if ($#ARGV==14) {
    $min=$ARGV[11];
    $max=$ARGV[12];
    $min_wave=$ARGV[13];
    $max_wave=$ARGV[14];
    $def=2;
}

open(FH,"<$pt");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	print "$x[$ns] $y[$ns]\n";
	$ns++;
    }
}
close(FH);

$pdl=rfits($unc_file);
($nx,$ny)=$pdl->dims;



open (OUT,">$outfile"); 
for ($i=0;$i<$ny;$i++) {
    $k=$i;
    print "Processing $i/$ny...";
    $call="img2spec.pl ".$unc_file." ".$i." spec_junk.txt";
    system($call);#
#c    print OUT "$i ";
    $call="compare_back_Z.pl spec_junk.txt ".$c_file." junk.out ".$mask_list." ".$vel_min." ".$vel_max." ".$delta_vel." ".$sigma." ".$wave_scale." ".$plot." ".$min." ".$max." ".$min_wave." ".$max_wave;
#    print "$call\n"; exit;
    system($call);
    open(T,"<COMPARE_BACK_Z.out");
    $line=<T>;
    chop($line);
    ($velocity,$chi_sq)=split(" ",$line);
    close(T);
    print OUT "$kc $x[$i] $y[$i] $velocity\n";
   # print "$i $x[$i] $y[$i] $velocity\n";
#		<stdin>;	      
}
close(OUT);

exit;
