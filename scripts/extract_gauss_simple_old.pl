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


if ($#ARGV<8) {
    print "USE: extract_gauss.pl RAW.fits Spectral_axis[0/1] TRACE.fits WIDTH OUTPUT.fits NX_MIN NX_MAX plot nplot [y_shift]\n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$trace_file=$ARGV[2];
$aperture=$ARGV[3];
$out_file=$ARGV[4];
$nx_min=$ARGV[5];
$nx_max=$ARGV[6];
$plot=$ARGV[7];
$nplot=$ARGV[8];
$y_shift=0;
if ($#ARGV==9) {
    $y_shift=$ARGV[9];
}
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

for ($i=0;$i<$nx_tr;$i++) {
    for ($j=0;$j<$ny_tr;$j++) {	
	$peak_y_max_new[$j][$i]=$peak_y_max_new[$j][$i]+$y_shift;
    }
}

print "Determining the parameters of the gaussians...\n";
for ($i=$nx_min;$i<$nx_max;$i++) {
#
# We create a cut for each peak
#
#$i=int($nx/2);
for ($k=0;$k<$ny_tr;$k++) {	
    $xk[$k]=$k;
    $range=int(2*$aperture);
    $j_min=$peak_y_max_new[$k][$i]-$range;
    $j_max=$peak_y_max_new[$k][$i]+$range;
	if ($j_min<0) {
	    $j_min=0;
	}
    if ($j_max>=$ny) {
	$j_max=$ny-1;
    }
    $min=1e12;
    $max=-1e12;
    my @a;
    my @xa;
    $na=0;
    for ($j=$j_min;$j<$j_max;$j++) {	
	$a[$na]=$in_array[$j][$i];
	if ($min>$a[$na]) {
	    $min=$a[$na];
	}
	if ($max<$a[$na]) {
	    $max=$a[$na];
	}
	$xa[$na]=$j;
	$na++;
    }
    $pdl_xa=pdl(@xa);
    $pdl_a=pdl(@a);
    ($pcen, $ppk, $pfwhm2, $pback, $perr, $pfit)=fitgauss1d($pdl_xa,$pdl_a);
    $cen[$k]=$pcen->slice(0)->sclr;
    $pk[$k]=$ppk->slice(0)->sclr;
    $fwhm2[$k]=$pfwhm2->slice(0)->sclr;
    $back[$k]=$pback->slice(0)->sclr;
    $err[$k]=$perr->slice(0)->sclr;

    if ($fwhm2[$k]>50) {
	$fwhm2[$k]=50;
    }

    $a_fwhm2[$k][$i]=$fwhm2[$k];
    $a_cen[$k][$i]=$cen[$k];
    $a_flux[$k][$i]=1.069*$pk[$k]*$fwhm2[$k];
    $a_pk[$k][$i]=$pk[$k];
    $a_back[$k][$i]=$back[$k];
    


#    print "($cen, $pk, $fwhm2, $back, $err, $fit)";
#    ($cen, $pk, $fwhm2, $back, $err, $fit) 
    if (($plot==1)&&($command ne "A")) {
	my @a_mod;
	for ($j=0;$j<$na;$j++) {
	    $a_mod[$j]=$pk[$k]*exp(-0.5*(($xa[$j]-$cen[$k])/($fwhm2[$k]/2.345))**2)+$back[$k]
	}


	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($j_min,$j_max,$min,$max,0,0);
	pglabel("Y-axis","Counts","$k/$ny_tr");
	pgpoint($na,\@xa,\@a,3);
	pgsci(8);
	pgline($na,\@xa,\@a_mod);
	pgclos();
#	$command="a";
	if ($command ne "a") {
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}
    }
}
#$command="P";
$med_fwhm2=median(@fwhm2);
$sig_fwhm2=sigma(@fwhm2);
if ($i==50*int($i/50)) {
    print "($i/$nx) FWHM=$med_fwhm2+-$sig_fwhm2\n";
}
#$npoly=3;
#$xk_pdl = pdl(@xk);
#$fwhm2_pdl = pdl(@fwhm2);
#($s_y,$coeff) = fitpoly1d $xk_pdl,$fwhm2_pdl,$npoly;
#for ($j=0;$j<$npoly;$j++) {
#    $c[$j]=$coeff->slice($j)->sclr;
#}
#for ($j=0;$j<$ny_tr;$j++) {
#    $fwhm2_out[$j]=0;
#    for ($kk=0;$kk<$npoly;$kk++) {
#	$fwhm2_out[$j]=$fwhm2_out[$j]+$c[$kk]*($j**$kk);
#    }
#}
@fwhm2_out=median_filter(6,\@fwhm2);
if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$ny_tr,$med_fwhm2-3*$sig_fwhm2,$med_fwhm2+3*$sig_fwhm2,0,0);
	pglabel("Spectral Peaks","FWHM","");
	pgpoint($ny_tr,\@xk,\@fwhm2,3);
	pgsci(2);
	pgline($ny_tr,\@xk,\@fwhm2_out);
	pgsci(8);
	pgline(2,[0,$ny_tr],[$med_fwhm2,$med_fwhm2]);
	pgclos();
	if ($command ne "a") {
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}
    }

#print "Press Enter"; <stdin>;
}

for ($i=0;$i<$nx_min;$i++) {
    for ($k=0;$k<$ny_tr;$k++) {	
	$a_fwhm2[$k][$i]=$a_fwhm2[$k][$nx_min];
	$a_cen[$k][$i]=$a_cen[$k][$nx_min];
	$a_flux[$k][$i]=$a_flux[$k][$nx_min];
	$a_back[$k][$i]=$a_back[$k][$nx_min];    
    }
}

for ($i=$nx_max;$i<$nx;$i++) {
    for ($k=0;$k<$ny_tr;$k++) {	
	$a_fwhm2[$k][$i]=$a_fwhm2[$k][$nx_max-1];
	$a_cen[$k][$i]=$a_cen[$k][$nx_max-1];
	$a_flux[$k][$i]=$a_flux[$k][$nx_max-1];
	$a_back[$k][$i]=$a_back[$k][$nx_max-1];    
    }
}

#$a_mfwhm2=median_smooth(@a_fwhm2,3,3);
print "DONE\n";
print "Writting fwhm.fits, cen.fits & bank.fits\n";
system("rm fwhm.fits");
#write_fits("fwhm.fits",[$nx,$ny_tr],2,\@a_mfwhm2);
write_fits("fwhm.fits",[$nx,$ny_tr],2,\@a_fwhm2);
system("rm cen.fits");
write_fits("cen.fits",[$nx,$ny_tr],2,\@a_cen);
system("rm back.fits");
write_fits("back.fits",[$nx,$ny_tr],2,\@a_back);
#system("rm flux.fits");
#write_fits("flux.fits",[$nx,$ny_tr],2,\@a_flux);
print "Writting the Aperture extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@a_flux);
exit;
