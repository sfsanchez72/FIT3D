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


if ($#ARGV<6) {
    print "USE: smooth_trace.pl TRACE.fits y_width Npoly OUT_TRACE.fits plot NX_min NX_max\n";
    exit;
}

$trace_file=$ARGV[0];
$y_width=$ARGV[1];
$npoly=$ARGV[2];
$out_file=$ARGV[3];
$plot=$ARGV[4];
$nx_min=$ARGV[5];
$nx_max=$ARGV[6];

print "Reading file $trace_file ";
$nax=read_naxes($trace_file);   
@naxis=@$nax;
@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx=$naxis[0];
$ny=$naxis[1];


for ($k=0;$k<$ny;$k++) {
    my @ax,@ay;
    $kk=0;
    for ($i=$nx_min;$i<$nx_max;$i++) {
	$ax[$kk]=$i;
	$ay[$kk]=$peak_y_max_new[$k][$i];
	$kk++;
    }
    $ax_pdl = pdl(@ax);
    $ay_pdl = pdl(@ay);
    ($s_y,$coeff) = fitpoly1d $ax_pdl,$ay_pdl,$npoly;
    for ($j=0;$j<$npoly;$j++) {
	$c[$j]=$coeff->slice($j)->sclr;
    }
    for ($i=0;$i<$nx;$i++) {
	$peak_y_max_out[$k][$i]=0;
	for ($j=0;$j<$npoly;$j++) {
	    $peak_y_max_out[$k][$i]=$peak_y_max_out[$k][$i]+$c[$j]*($i**$j);
	}
    }
}

print "DONE\n";
print "Writting the tracing solution $out_file...";
system("rm $out_file");
write_fits($out_file,[$nx,$ny],2,\@peak_y_max_out);
print "DONE\n";

if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    for ($k=0;$k<$ny;$k++) {
	pgenv(0,$nx,$peak_y_max_new[$k][int($nx/2)]-$y_width,$peak_y_max_new[$k][int($nx/2)]+$y_width,0,0);
	pglabel("X-axis","Y-Axis","$k/$ny");
	for ($i=0;$i<$nx;$i++) {
	    pgsci(4);
	    $x=$i;
	    $y=$peak_y_max_new[$k][$i];
	    pgpoint(1,[$x],[$y],5);
	    pgsci(1);
	    $x=$i;
	$y=$peak_y_max_out[$k][$i];
	    pgpoint(1,[$x],[$y],1);
	}
	pgsci(1);
	pgclose;
	if ($command ne "A") {
	    print "Press Enter"; 
	    $command=<stdin>; 
	    chop($command);
	} 
	if ($command eq "q") {
	    exit;
	}

    }
    pgend;
}

exit;
