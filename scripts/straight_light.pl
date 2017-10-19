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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<9) {
    print "USE: straight_light.pl RAW.fits Spectral_axis[0/1] TRACE.fits width Npoly staight_light.fits clean.fits NY_min NY_max plot\n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$trace_file=$ARGV[2];
$fwhm=$ARGV[3];
$npoly=$ARGV[4];
$out_file=$ARGV[5];
$clean_file=$ARGV[6];
$ny_min=$ARGV[7];
$ny_max=$ARGV[8];
$plot=$ARGV[9];
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
print "DONE\n";
print "Reading file $trace_file ";
$nax_tr=read_naxes($trace_file);   
@naxis_tr=@$nax_tr;
@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx_tr=$naxis_tr[0];
$ny_tr=$naxis_tr[1];

if ($nx_tr!=$nx) {
    print "The X-dimension of the trace file ($trace_file), $nx_tr\n";
    print "does not match that of the original image, $nx\n";
    exit;
}

print "Determining the Straight-light...\n";
open(CLEAN,">clean.txt");
for ($i=0;$i<$nx;$i++) {
    my @mask;
    for ($j=0;$j<$ny;$j++) {	
	$mask[$j]=1;
    }    
    for ($j=0;$j<$ny_tr;$j++) {
	$j_cen=int($peak_y_max_new[$j][$i]);
	$j_min=$j_cen-$fwhm;
	$j_max=$j_cen+$fwhm;
	for ($jj=$j_min;$jj<$j_max;$jj++) {
	    $mask[$jj]=0;
	}
    }
    $k=0;
    my @ax;
    my @ay;
    my @w;
    $min=1e12;
    $max=-1e12;
    $ny_med=int(0.5*($ny_min+$ny_max));
    for ($j=$ny_min;$j<$ny_max;$j++) {	
	if ($mask[$j]==1) {
	    $ax[$k]=$j;
	    $ii1=$i-2;
	    $ii2=$i+2;
	    if ($ii1<0) {
		$ii1=0;
		$ii2=5;
	    }
	    if ($ii2>$nx) {
		$ii1=$nx-5;
		$ii2=$nx;
	    }
	    $iiN=0;
#	    $ay[$k]=0;
#	    for ($ii=$ii1;$ii<$ii2;$ii++) {
#		$ay[$k]=$ay[$k]+$in_array[$j][$ii];
#		$iiN++;
#	    }
	$ay[$k]=$in_array[$j][$i];
#	    $ay[$k]=$ay[$k]/$iiN;
	    if ($min>$ay[$k]) {
		$min=$ay[$k];
	    }
	    if ($max<$ay[$k]) {
		$max=$ay[$k];
	    }
	    $w[$k]=(1+($j-$ny_med)**4);
	    print CLEAN "$k $i $j $ay[$k]\n";
	    $k++;
	}

    }

    @ay_tmp=median_filter(5,\@ay);
    @ay=@ay_tmp;


    

#
# We fit a polynomial function
#
    $ax_pdl = pdl(@ax);
    $ay_pdl = pdl(@ay);
    $w_pdl= pdl(@w);
    
    


    $YY=pdl([0..$ny]);
#    $YY=$YY+0.5;

    $pdl_sl_out=interpol($YY,$ax_pdl,$ay_pdl);


    ($s_y,$coeff) = fitpoly1d $YY,$pdl_sl_out,$npoly;
   
    for ($j=0;$j<$ny;$j++) {
	$sl_out[$j][$i]=$s_y->at($j);
    }

    @a_sl_out=list($pdl_sl_out);


    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$ny,$min,$max,0,0);
	pglabel("Y-axis","Counts","$i/$nx");
	pgpoint($k,\@ax,\@ay,3);
	for ($j=0;$j<$ny;$j++) {
	    $x=$j;
	    $y=$sl_out[$j][$i];
	    $yy=$a_sl_out[$j];
#	    pgsci(8);
#	    pgpoint(1,[$x],[$yy],16);
	    pgsci(2);
	    pgpoint(1,[$x],[$y],1);
	    pgsci(1);
	}


	print "Press Enter"; <stdin>;	
    }
    

}
close(CLEAN);
#print "COSA?";
print "DONE\n";
print "Smoothing...\n";
#print "Smoothing the initial values...\n";
$sum=0;
$nnx=31;
$nny=3;
$pdl_kernel=zeroes($nnx,$nny);
for ($j=0;$j<$nny;$j++) {
    for ($i=0;$i<$nnx;$i++) {
	$r=sqrt(($i-$nnx/2)**2+($j-$nny/2)**2);
	$val=exp(-0.5*($r/($nnx))**2);
	set($pdl_kernel,$i,$j,$val);
	$sum=$sum+$val;
#	print "$i $j $kernel[$i][$j]\n";
    }
}
#$pdl_kernel=pdl(@kernel);
$pdl_kernel=$pdl_kernel/$sum;
$pdl_out=pdl(@sl_out);
@dims=$pdl_out->dims;
#print "*=@dims\n";
#$smoothed=med2df($pdl_out,21,3,{Boundary => Reflect});
#$smoothed=med2df($pdl_fwhm2,21,21,{Boundary => Reflect});
$smoothed=conv2d($pdl_out,$pdl_kernel,{Boundary => Reflect});
#$smoothed2=med2df($smoothed,31,1,{Boundary => Reflect});
$smoothed2=med2df($smoothed,11,1,{Boundary => Reflect});
for ($j=0;$j<$ny;$j++) {	
    for ($i=0;$i<$nx;$i++) {
#	$msl_out[$j][$i]=$smoothed->at($i,$j);
	$msl_out[$j][$i]=$smoothed2->at($i,$j);
	$in_array[$j][$i]=$in_array[$j][$i]-$msl_out[$j][$i];
    }
}

#@msl_out=mean_smooth(\@sl_out,$nx,$ny,2,5);
print "DONE\n";
print "Writting the straight-light solution $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny],2,\@msl_out);
print "DONE\n";
print "Writting the cleaned data $clean_file\n";
system("rm $clean_file");
write_fits($clean_file,[$nx,$ny],2,\@in_array);
print "DONE\n";
exit;
