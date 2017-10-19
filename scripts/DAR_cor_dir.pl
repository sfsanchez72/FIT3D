#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<0) {
    print "DAR_cor_dir.pl filter[B,V]\n";
    exit;
}

$filter=$ARGV[0];
$n1=100;
$n2=1900;



open(DIR,"ls *.cube.fits |");
while($file=<DIR>) {
    chop($file);
    if (($file !~ "mask")&&($file !~ "mos")&&($file !~ "rad")&&($file !~ "std")) {
	print "$file ";
	$cube=$file;
	$cube_out=$file;
	$cube_out =~ s/cube/scube/;
	$V=$file;
	$V =~ s/cube/$filter/;
	$img=rfits($V);
	($nx,$ny)=$img->dims;
	$val_max=-1e12;
	$XC=38;
	$YC=42;
	for ($ii=35;$ii<45;$ii++) {
	    for ($jj=35;$jj<45;$jj++) {
		$val=$img->at($ii,$jj);
		if ($val>$val_max) {
		    $val_max=$val;
		    $XC=$ii;
		    $YC=$jj;
		}
	    }
	}
	print "$cube $XC $YC ";
	$call="rm -r ".$cube_out." ";
	system($call);
	$pdl=rfitshdr($cube);
	$nz=$pdl->{NAXIS3};
	print "NZ=$nz\n";
#	exit;

	if ($nz>1500) {
	    $call="DAR_det_cube.pl ".$cube."  ".$XC." ".$YC." 2 2 ".$cube_out." 200 1700 6 1 0 0 0 5";
	} else {
	    $call="DAR_det_cube.pl ".$cube."  ".$XC." ".$YC." 2 2 ".$cube_out." 70 1000 7 1 0 0 0 3";
	}
	print "$call\n"; exit;
	system($call);
	print "DONE\n";
	$name=$cube;
	$cut=".cube.fits";
	$name =~ s/$cut//;
	$XC=$XC+0.5;
	$YC=$YC+0.5;
	$call="radial_sum_cube.pl ".$name.".scube.fits 2.5 ".$XC." ".$YC." rad.rss.fits 2";
	system($call);
	$call="cp radial_sum_cube_0.ps rad_".$name."_5.ps";
	system($call);
	$call="img2spec.pl rad.rss.fits 0 ".$name.".spec_5.txt";
	system($call);
	$call="img2spec.pl rad.rss.fits 2 ".$name.".spec_10.txt";
	system($call);
	$call="cp radial_sum_cube_7.ps rad_".$name."_40.ps";
	system($call);
	$call="img2spec.pl rad.rss.fits 7 ".$name.".spec_40.txt";
	system($call);
	$call="cp rad.rss.fits rad_".$name.".rss.fits";
	system($call);
#	$call="kghostview radial_sum_cube_0.ps &";
#	system($call);


    }
}

