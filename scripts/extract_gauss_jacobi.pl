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


if ($#ARGV<4) {
    print "USE: extract_gauss_weight.pl RAW.fits Spectral_axis[0/1] TRACE.fits FWHM OUTPUT.fits [SHIFT] \n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
#$trace_file=$ARGV[2];
$cen_file=$ARGV[2];
$fwhm=$ARGV[3];
#$width=$ARGV[4];
#$shift=$ARGV[5];
#$shift_fwhm=$ARGV[6];
$out_file=$ARGV[4];
if ($#ARGV==5) {
    $y_shift=$ARGV[5];
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

#@a_fwhm2=read_img($fwhm_file);
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

$pdl_cen=pdl(@a_cen);
@a_cen=list($pdl_cen);
#$pdl_cen2=pdl(@a_cen);
#$pdl_cen2->wfits("mcen.fits");

print "DONE\n";
#
#
#
print "Gaussian Weighted, using the Jacobi Extraction\n";
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
	$a_flux[$jj][$i]=$in_array[int($a_cen[$jj][$i])][$i];
	for ($j=0;$j<$ny;$j++) {
	    $w[$jj][$j]=(exp(-0.5*(($j*1.0-$a_cen[$jj][$i])/($fwhm/2.345))**2))/(2.5066*($fwhm/2.345));
	}
    }

    for ($kk=1;$kk<$ny;$kk++) {
#	$guess=$a_flux[$kk][$i];
#	for ($ii=1;$ii<$ny_tr;$ii++) {
#	    $guess[$ii]=$a_flux[$ii][$i];
#	}
	for ($ii=1;$ii<$ny_tr;$ii++) {
	    do {	    		
		$xx=0;
		for ($jj=1;$jj<($ii-1);$jj++) {
		    $xx=$xx+$w[$jj][$kk]*$a_flux[$jj][$i];
		}
		for ($jj=($ii+1);$jj<$ny_tr;$jj++) {
		    $xx=$xx+$w[$jj][$kk]*$a_flux[$jj][$i];
		}
		$guess=($in_array[$kk][$i]-$a_flux[$ii][$i])/$w[$ii][$kk];		
		$delta=($guess-$a_flux[$ii][$i]);
		$a_flux[$ii][$i]=$guess;
	    } while ($delta>0.5);	    
	}
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

