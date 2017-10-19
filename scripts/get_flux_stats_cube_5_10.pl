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

#use POSIX;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<5) {
    print "USE:get_flux_stats.pl SPECTRUM.CUBE.FITS centra_WAVE width CRVAL CDELT DEVICE [AREA/arcsec^2]\n";    
    exit;
}

$spec_file=$ARGV[0];
$central_w=$ARGV[1];
$width=$ARGV[2];
$crval=$ARGV[3];
$cdelt=$ARGV[4];
$dev=$ARGV[5];
$area=1;
if ($#ARGV==6) {
    $area=$ARGV[6];
}

$pdl=rfits($spec_file);
($nz,$ny,$nx)=$pdl->dims();


$x_min=1e12;
$x_max=-1e12;
$y_min=1e12;
$y_max=-1e12;
$K=0;
$K3=0;
$K5=0;

$start_w=$central_w-0.5*$width;
$end_w=$central_w+0.5*$width;

$start_wb=$central_w-1*$width;
$end_wb=$central_w-0.5*$width;

$start_wr=$central_w+0.5*$width;
$end_wr=$central_w+1*$width;

$NT=$nz*$ny;

$J=0;
for ($ii=0;$ii<$nz;$ii++) {
    for ($jj=0;$jj<$ny;$jj++) {
	$k=0;
	$kc=0;
	my @f;
	my @fc;
	for ($i=0;$i<$nx;$i++) {
	    $w=$crval+$cdelt*$i;
	    if (($w>$start_w)&&($w<$end_w)) {
		$f[$k]=$pdl->at($ii,$jj,$i);
		$wave[$k]=$w;
		$k++;
	    }
	    if (($w>$start_wb)&&($w<$end_wb)) {
		$fc[$kc]=$pdl->at($ii,$jj,$i);
		if ($fc[$kc] ne "nan") {
		    $kc++;
		}
	    }
	    if (($w>$start_wr)&&($w<$end_wr)) {
		$fc[$kc]=$pdl->at($ii,$jj,$i);
#		print "KC = $w $fc[$kc] $kc\n";
		if ($fc[$kc] ne "nan") {
		    $kc++;
		}
	    }
	}
	$med[$J]=median(@f);
	$sig[$J]=sigma(@fc);
	$med_w=median(@wave);
	if ($sig[$J]>0) {
	    $SN[$J]=$med[$J]/$sig[$J];
	} else {
	    $SN[$J]=0;
	}
#	print "SN = $kc $med[$J] $sig[$J] $SN[$J]\n";

	if ($x_min>$med[$J]) {
	    $x_min=$med[$J];
	}
	if ($x_max<$med[$J]) {
	    $x_max=$med[$J];
	}
	if ($y_min>$SN[$J]) {
	    $y_min=$SN[$J];
	}
	if ($y_max<$SN[$J]) {
	    $y_max=$SN[$J];
	}
	if (($SN[$J]>=5)&&($SN[$J]<=8)) {
	    $F[$K]=$med[$J];
	    $K++;
	}

	if ($SN[$J]>=3) {
	    $K3++;
	}
	
	if ($SN[$J]>=5) {
	    $K5++;
	}
	if (($med[$J]>0)&&($sig[$J] ne "nan")) {
#	    print "$J $med[$J] $sig[$J] $SN[$J]\n";
#	    <stdin>;
	    $J++;
	} else {
#	    print "NONE\n";
	}

    }
}

$med_F=median(@F);
$sig_F=sigma(@F);

#  $f[$i]=((10**(-0.4*($mag[$i]-8.9)))*(1e-23)*(3e18))/($ww**2);

$mag1=8.9-2.5*log10(($med_F*1e-16*(($med_w)**2)/((1e-23)*(3e18)))/$area);

$mag=-21.109-2.5*log10($med_F*1e-16);


print "# N.points N.TOTAL Flux_limit sigma_limit med_wave magV_3sig N>3sigma N>5sigma\n";
print "# Flux in 1e-16 Erg s-1 cm-2 A-1) \n";
print "$K $NT $med_F $sig_F $med_w $mag $K3 $K5 $mag1\n";



pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
#pgenv($w[0],$w[$n-1],$min-5,$max+5,0,0);
#pgenv($x_min,$x_max,$y_min-1,$y_max+1,0,0);
pgenv(0,1.2*$x_max,0,1.5*$y_max+1,0,0);
pgsch(1.6);           # Set character height
pglabel("Median Flux ($start_w - $end_w)","S/N (simple estimation)","$spec_file mag_5s=$mag");
pgsch(2.2);           # Set character height
pgsci(1);
pgpoint($J,\@med,\@SN,16);
pgsci(2);
pgsls(2);
pgline(2,[0,1.2*$x_max],[3,3]);
pgsci(3);
pgsls(1);
pgline(2,[0,1.2*$x_max],[4,4]);
pgclose;
pgend;

exit;
