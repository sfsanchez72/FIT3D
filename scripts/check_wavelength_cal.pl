#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
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
    print "USE: check_wavelength.pl EXTRACTED.fits E_LINES.FILE CRVAL CDELT WIDTH_AA N.SEARCH PLOT\n";
    exit;
}

$input_cube=$ARGV[0];
$e_file=$ARGV[1];
$crval=$ARGV[2];
$cdelt=$ARGV[3];
$width=$ARGV[4];
$ns=$ARGV[5];
$plot=$ARGV[6];

$width=$width/$cdelt;
$k=0;
open(EFILE,"<$e_file");
open(OUT,">check_wavelength_cal.out");
print "#name input_wavelength median_wave mean_wave sigma_wave median_flux mean_flux sigma_flux\n";
print OUT "#name input_wavelength median_wave mean_wave sigma_wave median_flux mean_flux sigma_flux\n";
open(OUT2,">check_wavelength_cal.all");
open(OUT3,">ARC_2.disp.id");
while($eline=<EFILE>) {
    chop($eline);
    if ($eline !~ "#") {
	($wave,$name)=split(" ",$eline);
	$i_wave=($wave-$crval)/$cdelt;
#	print "$i_wave=($wave-$crval)/$cdelt\n"; exit;
	print OUT3 "$i_wave $wave\n";
	$call="check_eline.pl ".$input_cube." ".$i_wave." ".$width." ".$ns." ".$plot;
	system($call);
#	print "$call\n";	exit;
	open(IN,"<check_eline.out");
	$in=<IN>;
	chop($in);
	($med_w_peak,$mean_w_peak,$sig_w_peak,$med_flux_peak,$mean_flux_peak,$sig_flux_peak)=split(" ",$in);
	close(IN);
	$med_w_peak=$crval+$cdelt*$med_w_peak;
	$mean_w_peak=$crval+$cdelt*$mean_w_peak;
	$sig_w_peak=$cdelt*$sig_w_peak;
	print "$name $wave $med_w_peak $mean_w_peak $sig_w_peak $med_flux_peak $mean_flux_peak $sig_flux_peak\n";
	print OUT "$name $wave $med_w_peak $mean_w_peak $sig_w_peak $med_flux_peak $mean_flux_peak $sig_flux_peak\n";
	$shift[$k]=$wave-$med_w_peak;
	$s_shift[$k]=$sig_w_peak;
	$k++;


	open(IN,"<check_eline.all");
	while($in=<IN>) {
	    chop($in);
	    if ($in !~ "#") {
		$win=$in*$cdelt;
		print OUT2 "$win\n";
	    }
	}
	close(IN);
    }
}
close(OUT3);
close(OUT2);
close(EFILE);
$med_shift=median(@shift);
$mean_shift=mean(@shift);
$sig_shift=mean(@s_shift);
print "#NLINES .. median_shift mean_shift sigma_shift\n";
print "N.LINES $k $med_shift $mean_shift $sig_shift\n";
print OUT "#NLINES .. median_shift mean_shift sigma_shift\n";
print OUT "N.LINES $k $med_shift $mean_shift $sig_shift\n";
close(OUT);



exit;

exit; 
