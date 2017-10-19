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


if ($#ARGV<10) {
    print "USE: calib_det.pl UNCALIBRATED_SPEC.txt CALIBRATED_SPEC.txt ratio.txt FACTOR EXTINC_AV AIRMASS SMOOTH_BOX1 SMOOTH2_BOX plot min_wave max_wave\n";
    exit;
}
$unc_file=$ARGV[0];
$c_file=$ARGV[1];
$out_file=$ARGV[2];
$factor=$ARGV[3];
$ext_AV=$ARGV[4];
$airmass=$ARGV[5];
$box=$ARGV[6];
$smooth=$ARGV[7];
$plot=$ARGV[8];
$min_wave=$ARGV[9];
$max_wave=$ARGV[10];

$n_unc=0;
open(FH,"<$unc_file");
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
	$index_unc[$n_unc]=$data[0];
	$wave_unc[$n_unc]=$data[1];
	$Krc=0.085*($wave_unc[$n_unc]/5450)**(-4);
#	$ext=1.1*$Krc+(0.8*$ext_AV*$airmass-1.1*0.085)*($wave_unc[$n_unc]/5450)**(-0.8);
	$ext=1.1*$Krc+(0.8*$ext_AV-1.1*0.085)*($wave_unc[$n_unc]/5450)**(-0.8);
#	$ext=$ext_AV*($wave_unc[$n_unc]/5450)**(-1.2);
	$Fext=10**(-0.4*$ext*$airmass);
	$flux_unc_ini[$n_unc]=$data[2]*$factor/($Fext);
	$n_unc++;
    }
}
close(FH);
$dpix_unc=$wave_unc[1]-$wave_unc[0];
print "Wavelength Range=$wave_unc[0] $wave_unc[$n_unc-1]\n";
if ($box>2) {
    @flux_unc=median_filter($box,\@flux_unc_ini);
} else {
#    if ($box<2) {
#	@flux_unc=median_filter(abs($box),\@flux_unc_ini);
#    } else {
	for ($i=0;$i<($#flux_unc_ini+1);$i++) {
	    $flux_unc[$i]=$flux_unc_ini[$i];
	}
#    }
}
$n_c=0;
open(FH,"<$c_file");
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
	$wave_c[$n_c]=$data[0];
	$flux_c_ini[$n_c]=$data[1];
	$n_c++;
    }
}
$dpix_c=$wave_c[1]-$wave_c[0];
close(FH);
if ($box>2) {
    @flux_c=median_filter($box,\@flux_c_ini);
} else {
    for ($i=0;$i<($#flux_c_ini+1);$i++) {
	$flux_c[$i]=$flux_c_ini[$i];
    }
}

print "$n_unc $n_c\n";
#
# We interpolate
#
$out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), pdl(@flux_c));
#print "$out_spec_pdl\n";
$min=1e12;
$max=-1e12;
$min_rat=1e30;
$max_rat=-1e30;
for ($j=0;$j<$n_unc;$j++) {
    $out_spec[$j]=$out_spec_pdl->at($j);		
# Error in the flux calibration
#
#    $rat1=($dpix_c/$dpix_unc);    


#    print "$out_spec[$j] $flux_unc[$j]\n";
    if ($flux_unc[$j]>0) {
	$rat2=($out_spec[$j]/$flux_unc[$j]);
    } else {
	$rat2=1e20;
    }
    $rat1=1;
    $ratio[$j]=$rat2/$rat1;
#    print "$j $ratio[$j]=$out_spec[$j]/$flux_unc[$j]\n";
#    if (($j>0.05*$n_unc)&&($j<0.9*$n_unc)) {
	if ($max<$flux_unc[$j]) {
	    $max=$flux_unc[$j];
	}
	if ($min>$flux_unc[$j]) {
	    $min=$flux_unc[$j];
	}
#    }

#    print "$j $ratio[$j] $min_rat $max_rat\n";
#    print "$j $wave_unc[$j] $flux_unc[$j] $out_spec[$j] $ratio[$j]\n";
}

#$max_rat=100;
$median=median(@ratio);
$i_max_wave=int(($max_wave-$wave_unc[0])/($wave_unc[1]-$wave_unc[0]));
$i_min_wave=int(($min_wave-$wave_unc[0])/($wave_unc[1]-$wave_unc[0]));
for ($j=0;$j<$n_unc;$j++) {
    $flux_unc_plot[$j]=$median*$flux_unc[$j];
    if ($wave_unc[$j]>$max_wave) {
	$ratio[$j]=$ratio[$i_max_wave];
    }
    if ($wave_unc[$j]<$min_wave) {
	$ratio[$j]=$ratio[$i_min_wave];
    }
	if ($max_rat<$ratio[$j]) {
	    $max_rat=$ratio[$j];
	}
	if ($min_rat>$ratio[$j]) {
	    $min_rat=$ratio[$j];
	}

}
#$min_rat=0.1*$median;
#$max_rat=0.8*$max_rat;


for ($j=0;$j<$n_c;$j++) {
    if (($wave_c[$j]>$wave_unc[0])&&($wave_c[$j]<$wave_unc[$n_unc-1])) {
	if ($max<$flux_c[$j]) {
	    $max=$flux_c[$j];
	}
	if ($min>$flux_c[$j]) {
	    $min=$flux_c[$j];
	}
    }
}


#$mratio=med2df(pdl(@ratio),31,1,{Boundary => Reflect});

#$npoly=2;
#@x_ratio=median_filter($smooth,\@wave_unc);
#@y_ratio=median_filter($smooth,\@ratio);

if ($smooth>2) {
    @x_med=median_box($smooth,\@wave_unc);
    @y_med=median_box($smooth,\@ratio);
} else {
    for ($i=0;$i<($#wave_unc+1);$i++) {
	$x_med[$i]=$wave_unc[$i];
	$y_med[$i]=$ratio[$i];
#	print "$i $x_med[$i]=$wave_unc[$i] $y_med[$i]=$ratio[$i]\n";
    }
}



for ($i=0;$i<($#wave_unc+1);$i++) {
    if (($wave_unc[$i]>6850)&&($wave_unc[$i]<8000)) {
#	$x_med[$i]=$wave_unc[$i];
#	$y_med[$i]=$ratio[$i];
    }
#	print "$i $x_med[$i]=$wave_unc[$i] $y_med[$i]=$ratio[$i]\n";
}

#************NEW

my $spline=new Math::Spline(\@x_med,\@y_med);
for ($i=0;$i<$n_unc;$i++) {
    $x_ratio[$i]=$wave_unc[$i];
    $y_ratio[$i]=$spline->evaluate($x_ratio[$i]);
}


#($pdl_s_ratio,$coeff) = fitpoly1d(pdl(@x_ratio),pdl(@y_ratio),$npoly);
#for ($j=0;$j<$#x_ratio+1;$j++) {
#    $s_ratio[$j]=$pdl_s_ratio->at($j);
#}


if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsubp(1,2);
    pgenv($wave_unc[0],$wave_unc[$n_unc-1],$min,$max,0,0);
    pgsch(1.2);           # Set character height
    pglabel("Wavelength","Counts","");
#    pgsch(0.9);           # Set character height
    pgsci(1);
    pgpoint($n_c,\@wave_c,\@flux_c,2);    
    pgsci(4);
    pgpoint($n_unc,\@wave_unc,\@flux_unc_ini,5);    
    pgsci(3);
#    pgline($n_unc,\@wave_unc,\@flux_unc);    
    pgline($n_unc,\@wave_unc,\@flux_unc_plot);    
    pgsci(2);
    pgline($n_unc,\@wave_unc,\@out_spec);    
    pgsci(1);
#    pgsci(3);
#    pgpoint($n_unc,\@wave_unc,\@ratio,5);    
     pgenv($wave_unc[0],$wave_unc[$n_unc-1],$min_rat,$max_rat,0,0);
#    pgsch(1.6);           # Set character height
    pglabel("Wavelength","Ratio","");
    pgsch(0.9);           # Set character height
    pgsci(1);
    pgpoint($n_unc,\@wave_unc,\@ratio,2);       
    pgsci(2);
#    pgline($n_unc,\@wave_unc,\@s_ratio);       
    pgline($#x_ratio,\@x_ratio,\@y_ratio);       
#    pgsci(5);
#    pgline($n_unc,\@wave_unc,\@s_ratio);       
#    pgline($#x_ratio,\@x_ratio,\@s_ratio);       

    pgsci(1);
    pgclose;
    pgend;
    print "Press Enter"; <stdin>;
}

open(FH,">$out_file");
for ($j=0;$j<$n_unc;$j++) {
    print FH "$j $wave_unc[$j] $y_ratio[$j]\n";
}
close(FH);

exit;

