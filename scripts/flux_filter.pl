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
if ($#ARGV<2) {
    print "USE: flux_filter.pl filter.dat spectrum redshift \n";
    exit;
}

$C=299792.458;

$filter=$ARGV[0];
#$mag_zero=$ARGV[1];
$input=$ARGV[1];
$redshift=$ARGV[2];
$vel_mod=$redshift*$C;

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
    $n++;
}
close(FH);

#<stdin>;

$nf=0;
open(FH,"<$filter");
while ($line=<FH>) {
    chop($line);
    if ($line != "") {
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
	if ($flux_f[$nf]<0) {
	    $flux_f[$nf]=0;
	}
#    $flux_f[$nf]=$data[2];
#    $wave_f[$nf]=$data[1];

#    print "$nf $flux_f[$nf] $wave_f[$nf] $waver_f[$nf]\n";
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
for ($i=0;$i<$n;$i++) {
    $val=$pdl_flux->at($i);
    $val=$val*$wave[$i];
    $new_val=$flux[$i]*$val;
    if (($val>0)&&($new_val>0)) {
	$Flux=$Flux+$new_val;
	$w=$w+$val;
    }
    $val=$pdl_flux_r->at($i);
     $val=$val*$wave[$i];
    $new_val=$flux[$i]*$val;
    if (($val>0)&&($new_val>0)) {
	$Flux_r=$Flux_r+$new_val;
	$w_r=$w_r+$val;
    }
#    print "$i $wave[$i] $flux[$i] $new_val $Flux $w $Flux_r $w_r\n";
}

#$Flux=($Flux/$w)*($wave[1]-$wave[0]);
if ($w>0) {
    $Flux=($Flux/$w);#*($wave[1]-$wave[0]);
}
if ($w_r>0) {
#    $Flux_r=$Flux_r/$w_r*($wave[1]-$wave[0]);
    $Flux_r=$Flux_r/$w_r;#*($wave[1]-$wave[0]);
}

open(FH,">flux_filter.out");
print FH "$Flux_r\n";
close(FH);
print "$Flux_r\n";

exit;
