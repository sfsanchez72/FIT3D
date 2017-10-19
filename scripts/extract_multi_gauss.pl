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

#
# We create a cut
#
    my @cut;
    my @xcut;
    for ($j=0;$j<$ny;$j++) {
	$cut[$j]=$in_array[$j][$i];
	$xcut[$j]=$j;
	if ($cut[$j]!=0) {
	    $ecut[$j]=1/$cut[$j]**2;
	} else {
	    $ecut[$j]=1;
	}
    }

#
# We create the models
#
    $A=0.5*($FWHM/2.345)**2;
    $NORM=sqrt(3.1416/$A);

    $n_mod=$ny_tr+1;
    $pdl_mod=zeroes($ny,$n_mod);
    for ($k=0;$k<$ny_tr;$k++) {	
	for ($j=0;$j<$ny;$j++) {
	    $a_mod=exp(-0.5*(($j-$peak_y_max_new[$k][$i])/($FWHM/2.345))**2)/$NORM;
	    set($pdl_mod,$j,$k,$a_mod);
	}
    }
#    print "PP\n";
    for ($j=0;$j<$ny;$j++) {
	set($pdl_mod,$j,$ny_tr,1);
    }
    $pdl_flux=pdl(@cut);
    $pdl_error=pdl(@ecut);	
#	$pdl_mod=pdl(@a_mod);
    ($yfit, $coeffs) = linfit1d $pdl_flux,$pdl_mod;;#,$pdl_error;        
    @out_fit=list($yfit);
#	print "$flux[$k][$i] $lcoeffs[0] $lcoeffs[1]\n";
    for ($k=0;$k<$ny_tr;$k++) {	
	$FLUX[$k][$i]=$coeffs->at($k);
    }
    $back=$coeffs->at($ny_tr);
    ($min,$max)=minmax(@cut);
    $nplot=4;
    if ($plot==1) {	
	pgbegin(0,"/xs",1,1);
	pgsubp(1,$nplot);
	for ($I=0;$I<$nplot;$I++) {
	    pgenv(($ny/$nplot)*$I,($ny/$nplot)*($I+1),$min,$max,0,0);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
#	    pgenv(0,$ny,$min,$max,0,0);
	    pglabel("Y-axis","Counts","$k/$ny_tr");
	    pgpoint($ny,\@xcut,\@cut,3);
	    pgsci(8);
	    pgline($ny,\@xcut,\@out_fit);
	    pgsci(1);

	}
	pgclos();
#	 pgend;
	print "BACK= $back, Press Enter"; 
	$command=<stdin>;
	    chop($command);
    }
   
}
print "DONE\n";
print "DONE\n";
print "Writting the gaussian extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@FLUX);
print "DONE\n";
exit;
