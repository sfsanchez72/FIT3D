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
#use PDL::Fit::Linfit;
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<6) {
    print "USE: extract_gauss_weight.pl RAW.fits Spectral_axis[0/1] CEN.fits FWHM.fits ADDING_WIDTH ALLOWED_SHIFT_CEN OUTPUT.fits [SHIFT] \n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
#$trace_file=$ARGV[2];
$cen_file=$ARGV[2];
$fwhm_file=$ARGV[3];
$width=$ARGV[4];
$shift=$ARGV[5];
#$shift_fwhm=$ARGV[6];
$out_file=$ARGV[6];
if ($#ARGV==7) {
    $y_shift=$ARGV[7];
}



#$nx_min=$ARGV[6];
#$nx_max=$ARGV[7];
#$plot=$ARGV[6];
#$nplot=$ARGV[7];
$command="P";
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
print "Reading the Centroid an FWHM files\n";

@a_fwhm2=read_img($fwhm_file);
@a_cen=read_img($cen_file);
$nax_tr=read_naxes($fwhm_file);   
@naxis_tr=@$nax_tr;
#@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx_tr=$naxis_tr[0];
$ny_tr=$naxis_tr[1];

for ($i=0;$i<$nx_tr;$i++) {
    for ($j=0;$j<$ny_tr;$j++) {	
	$a_cen[$j][$i]=$a_cen[$j][$i]+$y_shift;
    }
}

if ($nx_tr!=$nx) {
    print "The X-dimension of the trace file ($trace_file), $nx_tr\n";
    print "does not match that of the original image, $nx\n";
    exit;
}

print "Smoothing the initial values...\n";
$sum=0;
for ($j=0;$j<11;$j++) {
    for ($i=0;$i<11;$i++) {
	$r=sqrt(($i-4)**2+($j-4)**2);
	$kernel[$i][$j]=exp(-0.5*($r/3)**2);
	$sum=$sum+$kernel[$i][$j];
#	print "$i $j $kernel[$i][$j]\n";
    }
}
$pdl_kernel=pdl(@kernel);
$pdl_kernel=$pdl_kernel/$sum;

#$pdl_kernel=pdl(@kernel);

$pdl_fwhm2=pdl(@a_fwhm2);
#$smoothed=med2df($pdl_fwhm2,5,5,{Boundary => Default});
for ($i=0;$i<$nx;$i++) {
    my @ar;
    for ($j=0;$j<$ny_tr;$j++) {	
	$ar[$j]=$a_fwhm2[$j][$i];
    }
    $mean=median(@ar);
    $sigma=sigma(@ar);
    for ($j=0;$j<$ny_tr;$j++) {	
	if (abs($ar[$j]-$mean)>$sigma) {
	    $a_fwhm2[$j][$i]=$mean;
	}
    }
}

$smoothed=med2df($pdl_fwhm2,21,21,{Boundary => Reflect});
$smoothed2=conv2d($smoothed,$pdl_kernel,{Boundary => Reflect});
$smoothed3=med2df($smoothed2,5,5,{Boundary => Reflect});
for ($j=0;$j<$ny_tr;$j++) {	
    for ($i=0;$i<$nx;$i++) {
	$a_fwhm2[$j][$i]=$smoothed3->at($i,$j);
#	print "$i,$j $a_fwhm2[$j][$i]\n";
	}
}
$smoothed->wfits("mfwhm.fits");
$pdl_cen=pdl(@a_cen);
$smoothed=med2df($pdl_cen,3,3,{Boundary => Reflect});
#$smoothed2=conv2d($smoothed,$pdl_kernel,{Boundary => Reflect});
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny_tr;$j++) {	
	$val=$smoothed->at($i,$j);
	$delta=abs($val-$a_cen[$j][$i]);
#	print "$i,$j $delta $shift\n";
	if ($delta>$shift) {
#	    if (($j<2)||($j>$ny_tr)) {
		$a_cen[$j][$i]=$smoothed->at($i,$j);
#	    } else {
#		$a_cen[$j][$i]=0.5*($a_cen[$j-2][$i]+$a_cen[$j+2][$i]);
#	    }
	}
    }
}

#for ($i=0;$i<$nx;$i++) {
#    for ($j=1;$j<$ny_tr-1;$j++) {	
#	if (!(($a_cen[$j-1][$i]<$a_cen[$j][$i])&&($a_cen[$j][$i]<$a_cen[$j+1][$i]))) {
#	    $a_cen[$j][$i]=0.5*($a_cen[$j-1][$i]+$a_cen[$j+1][$i]);
#	}
#    }
#}
$pdl_cen2=pdl(@a_cen);
$pdl_cen2->wfits("mcen.fits");



print "DONE\n";
#
#
#
print "Gaussian Weighted Extraction\n";
#$command="P";
#$plot=1;
for ($i=0;$i<$nx;$i++) {
#for ($i=$nx_min;$i<$nx_max;$i++) {
    $min=1e12;
    $max=-1e12;
    my @a;
    my @sig;
    my @xa;
    my @cen;
    my @fwhm;
    my @back;
    $pdl_val=zeroes($ny,$ny_tr);
    for ($jj=0;$jj<$ny_tr;$jj++) {
	$j_min=int($a_cen[$jj][$i])-$width;
	$j_max=int($a_cen[$jj][$i])+$width;
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($j_max>$ny) {
	    $j_max=$ny;
	}
	$a_flux[$jj][$i]=0;
	$sum_w=0;
	for ($j=$j_min;$j<$j_max;$j++) {
	    $w=exp(-0.5*(($j*1.0-$a_cen[$jj][$i])/($a_fwhm2[$jj][$i]/2.345))**2);
	    $a_flux[$jj][$i]=$a_flux[$jj][$i]+$w*$in_array[$j][$i];
	    $sum_w=$sum_w+$w;
#	    print "$a_flux[$jj][$i]+$w*$in_array[$j][$i]\n";
	}
	$a_flux[$jj][$i]=$a_flux[$jj][$i]/$sum_w;
    }
    if ($i==50*int($i/50)) {
	print "$i/$nx\n";
    }
}
print "\n";

print "Writting the gaussian-weighted extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@a_flux);




print "DONE\n";
exit;


sub my_gauss1d {
    my $x=@_[0];
    my $cen=@_[1];
    my $fwhm=@_[2];
    my $valgauss=exp(-0.5*(($x-$cen)/($fwhm/2.345))**2);
    return $valgauss;
}

