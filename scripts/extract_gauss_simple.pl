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
use  PDL::Fit::Linfit;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<6) {
    print "USE: extract_aper.pl RAW.fits Spectral_axis[0/1] TRACE.fits Aperture OUTPUT.fits WIDTH PLOT [SHIFT] \n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$trace_file=$ARGV[2];
$aperture=$ARGV[3];
$out_file=$ARGV[4];
$FWHM=$ARGV[5];
$plot=$ARGV[6];
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
    if ($i==50*int($i/50)) {
	print "$i/$nx\n";
    }
    for ($k=0;$k<$ny_tr;$k++) {	
#	$range=int(2*$aperture);
	$range=(1.2*$aperture);
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
	my @cut;
	my @xcut;
	my $ncut=0;
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
	    $cut[$ncut]=$in_array[$j][$i];
	    if ($cut[$ncut]!=0) {
		$ecut[$ncut]=1/$cut[$ncut]**2;
	    } else {
		$ecut[$ncut]=1;
	    }
	    $xcut[$ncut]=$j;
	    $ncut++;
	}

#
# We prepare the data
#

	#
	# We removed the adjacent spaxels
	# 
	$A=0.5*($FWHM/2.345)**2;
	$NORM=sqrt(3.1416/$A)/1.5;

	if ($k>0) {
	    for ($kk=0;$kk<$ncut;$kk++) {
		$G=exp(-0.5*(($xcut[$kk]-$peak_y_max_new[$k-1][$i])/($FWHM/2.345))**2)/$NORM;
		$cut[$kk]=$cut[$kk]-$flux[$k-1][$i]*$G;
	    }
	}

	if ($k<$ny_tr-1) {
	    for ($kk=0;$kk<$ncut;$kk++) {
		$G=exp(-0.5*(($xcut[$kk]-$peak_y_max_new[$k+1][$i])/($FWHM/2.345))**2)/$NORM;
		$cut[$kk]=$cut[$kk]-$flux[$k+1][$i]*$G;
	    }
	}

	$pdl_flux=pdl(@cut);
	$pdl_error=pdl(@ecut);
	$pdl_mod=zeroes($ncut,2);
	for ($kk=0;$kk<$ncut;$kk++) {
	    $a_mod[$kk]=exp(-0.5*(($xcut[$kk]-$peak_y_max_new[$k][$i])/($FWHM/2.345))**2)/$NORM;
	    set($pdl_mod,$kk,0,$a_mod[$kk]);
	    set($pdl_mod,$kk,1,1);
	}
	
	
#	$pdl_mod=pdl(@a_mod);
	($yfit, $coeffs) = linfit1d $pdl_flux,$pdl_mod;#,$pdl_error;
	@out_fit=list($yfit);
	@lcoeffs=list($coeffs);
#	print "$flux[$k][$i] $lcoeffs[0] $lcoeffs[1]\n";
	$FLUX[$k][$i]=$lcoeffs[0];
	($min,$max)=minmax(@cut);
	if ($plot==1) {

	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv($j_min,$j_max,$min,$max,0,0);
	    pglabel("Y-axis","Counts","$k/$ny_tr");
	    pgpoint($ncut,\@xcut,\@cut,3);
	    pgsci(8);
	    pgline($ncut,\@xcut,\@out_fit);
	    pgclos();
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}

    }
}
print "DONE\n";
print "DONE\n";
print "Writting the gaussian extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@FLUX);
print "DONE\n";
exit;
