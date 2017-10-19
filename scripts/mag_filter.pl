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

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{'PGPLOT_ENVOPT'}="V";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
if ($#ARGV<3) {
    print "USE: mag_filter.pl filter.dat zero_point spectrum redshift [WMIN WMAX]\n";
    exit;
}

$C=299792.458;

$filter=$ARGV[0];
$mag_zero=$ARGV[1];
$input=$ARGV[2];
$redshift=$ARGV[3];
$vel_mod=$redshift*$C;

$wmin=0;
$wmax=1e14;
if ($#ARGV==5) {
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
}

$n=0;
open(FH,"<$input");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {

	$flux[$n]=$data[2];
	$wave[$n]=$data[1];

    } else {
	$flux[$n]=$data[1];
	$wave[$n]=$data[0];
    }

#    $wave[$n]=$data[1];
#    $flux[$n]=$data[2];
#    print "$n $wave[$n] $flux[$n]\n";
    if ($flux[$n]>0) {
	if (($wave[$n]>$wmin)&&($wave[$n]<$wmax)) {
	    $n++;
	}
    }
}
close(FH);

#<stdin>;

$nf=0;
open(FH,"<$filter");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
	$flux_f[$nf]=$data[2];
	$wave_f[$nf]=$data[1];
	$waver_f[$nf]=$data[1]/(1+$vel_mod/$C);
    } else {
	$flux_f[$nf]=$data[1];
	$wave_f[$nf]=$data[0];
	$waver_f[$nf]=$data[0]/(1+$vel_mod/$C);
    }
#    $flux_f[$nf]=$data[2];
#    $wave_f[$nf]=$data[1];

#    print "$nf $flux_f[$nf] $wave_f[$nf] $waver_f[$nf]\n";
    if ($flux_f[$nf]>0) {
	$nf++;
    }
}
close(FH);

#<stdin>;
#print "PASO\n";
my $pdl_flux = interpol(pdl(@wave), pdl(@wave_f), pdl(@flux_f));
my $pdl_flux_r = interpol(pdl(@wave), pdl(@waver_f), pdl(@flux_f));
#print "PASO 2\n";

$Flux=0;
$Flux_r=0;
$w=0;
$w_r=0;
$sum++;
for ($i=0;$i<$n;$i++) {
    $val=$pdl_flux->at($i);
    $val_F=$pdl_flux->at($i);
    $val=$val*$wave[$i];
    $new_val=$flux[$i]*$val;
#    print "$i $wave[$i] $val $new_val $flux[$i]\n";
    if ($val>0) {
	$Flux=$Flux+$new_val;
	$w=$w+$val;
#	$W=$W+$wave[$i]
	$sum=$sum+$val_F;
    }
    $val=$pdl_flux_r->at($i);
    $val=$val*$wave[$i];
    $new_val=$flux[$i]*$val;
    if ($val>0) {
	$Flux_r=$Flux_r+$new_val;
	$w_r=$w_r+$val;
    }
#    print "$i $wave[$i] $flux[$i] $new_val $Flux $w $Flux_r $w_r\n";
}
$w_eff=$w/$sum;

$Flux=($Flux/$w);#*($wave[1]-$wave[0]);#*($wave_f[1]-$wave_f[0]);
#$Flux=($Flux)#*($wave[1]-$wave[0]);
$Flux_r=$Flux_r/$w_r;#*($wave[1]-$wave[0]);
#$Flux_norm=$Flux/$w_eff;

$mag=$mag_zero-2.5*log10($Flux);
$mag_r=$mag_zero-2.5*log10($Flux_r);
$mag_AB=8.9-2.5*log10(($Flux*1e-16*(($w_eff)**2)/((1e-23)*(3e18))));


print "$input $mag $mag_r $w_eff $mag_AB\n";

exit;
