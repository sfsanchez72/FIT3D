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
    print "USE: extract_aper_CT.pl RAW.fits Spectral_axis[0/1] TRACE.fits Aperture OUTPUT.fits FWHM delta [SHIFT]\n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$trace_file=$ARGV[2];
$aperture=$ARGV[3];
$out_file=$ARGV[4];
$fwhm=$ARGV[5];
$delta=$ARGV[6];
$y_shift=0;
if ($#ARGV==7) {
    $y_shift=$ARGV[7];
}



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
for ($i=0;$i<$nx_tr;$i++) {
    for ($j=0;$j<$ny_tr;$j++) {	
	$peak_y_max_new[$j][$i]=$peak_y_max_new[$j][$i]+$y_shift;
    }
}


if ($nx_tr!=$nx) {
    print "The X-dimension of the trace file ($trace_file), $nx_tr\n";
    print "does not match that of the original image, $nx\n";
    exit;
}

$a=0.5*$aperture;
print "Extracting...\n";
for ($i=0;$i<$nx;$i++) {
#    if ($i==50*int($i/50)) {
	print "$i/$nx\n";
#    }
    for ($k=0;$k<$ny_tr;$k++) {	
	$range=int(1.5*$aperture);
	$j_min=int($peak_y_max_new[$k][$i]-$range);
	$j_max=int($peak_y_max_new[$k][$i]+$range);
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($j_max>=$ny) {
	    $j_max=$ny-1;
	}
	$flux[$k][$i]=0;
	$sum_w=0;
	$q=$peak_y_max_new[$k][$i];
	for ($j=$j_min;$j<$j_max;$j++) {
	    $w=0;
	    $qa1=$q-$a;
	    $qa2=$q+$a;
	    $p1=$j-0.5;
	    $p2=$j+0.5;
	    if ((($q-$a)<=($j-0.5))&&(($q+$a)>=($j+0.5))) {
		$w=1;
	    }
	    if ((($q-$a)<=($j-0.5))&&(($q+$a)<($j+0.5))) {
		$w=$q+$a-($j-0.5);
	    }
	    if ((($q-$a)>($j-0.5))&&(($q+$a)>=($j+0.5))) {
		$w=$j+0.5-($q-$a);
	    }
	    if ($w<0) {
		$w=0;
	    }
	    $sum_w=$sum_w+$w;
	    $flux[$k][$i]=$flux[$k][$i]+$w*$in_array[$j][$i];
	}
	for ($j=0;$j<$ny;$j++) {
	    $we[$k][$j]=(exp(-0.5*(($j*1.0-$peaks_y_max_new[$k][$i])/($fwhm/2.345))**2))/(2.5066*($fwhm/2.345));
	}
	$conver[$k]=0;
	$n_conver=0;
    }

#    print "It";
 
#
#Iteration
#
    $niter=0;
    while (($n_conver<$ny_tr)&&($niter<5*$ny_tr)) {	
	for ($k=0;$k<$ny_tr;$k++) {	
	    if ($conver[$k]==0) {
		$range=int(1.5*$aperture);
		$j_min=int($peak_y_max_new[$k][$i]-$range);
		$j_max=int($peak_y_max_new[$k][$i]+$range);
		if ($j_min<0) {
		    $j_min=0;
		}
		if ($j_max>=$ny) {
		    $j_max=$ny-1;
		}
		$flux_new[$k][$i]=0;
		$sum_w=0;
		$q=$peak_y_max_new[$k][$i];
		for ($j=$j_min;$j<$j_max;$j++) {
		    $w=0;
		    $qa1=$q-$a;
		    $qa2=$q+$a;
		    $p1=$j-0.5;
		    $p2=$j+0.5;
		    if ((($q-$a)<=($j-0.5))&&(($q+$a)>=($j+0.5))) {
			$w=1;
		    }
		    if ((($q-$a)<=($j-0.5))&&(($q+$a)<($j+0.5))) {
			$w=$q+$a-($j-0.5);
		    }
		    if ((($q-$a)>($j-0.5))&&(($q+$a)>=($j+0.5))) {
			$w=$j+0.5-($q-$a);
		    }
		    if ($w<0) {
			$w=0;
		    }
		    $sum_w=$sum_w+$w;		
		    $flux_new[$k][$i]=$flux_new[$k][$i]+$w*$in_array[$j][$i];
		    if ($k>0) {
			$flux_new[$k][$i]=$flux_new[$k][$i]-$we[$k-1][$j]*$flux[$k-1][$i];
		    }
		    if ($k<($ny_tr-1)) {
		    $flux_new[$k][$i]=$flux_new[$k][$i]-$we[$k+1][$j]*$flux[$k+1][$i];
		}		    
		
		}
	    }
	    $niter++;
#	    print "$k,$i $flux_new[$k][$i]-$flux[$k][$i])<$delta\n";
	    if (abs($flux_new[$k][$i]-$flux[$k][$i])<$delta) {
		$conver[$k]=1;
		$n_conver++;
	    } 	
	}	
#       @flux_new=@flux;
	@flux=@flux_new;

#	print "$n_conver, $niter\n";
    }
#    print 
}




print "DONE\n";
print "DONE\n";
print "Writting the Aperture extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@flux);
print "DONE\n";
exit;
