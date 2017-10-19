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
    print "USE:get_flux_width_cube.pl SPECTRUM.FITS ELINES.LIST WIDTH CRVAL CDELT DEV\n";    
    exit;
}

$spec_file=$ARGV[0];
$e_file=$ARGV[1];
$width=$ARGV[2];
$crval=$ARGV[3];
$cdelt=$ARGV[4];
$dev=$ARGV[5];

$pdl_cube=rfits($spec_file);
($NX,$NY,$NZ)=$pdl_cube->dims();

$nz=int($NZ/2);
$K=0;
for ($i=0;$i<$NX;$i++) {
    for ($j=0;$j<$NY;$j++) {
	$val= $pdl_cube->at($i,$j,$nz);
	if ($val>0) {
	    $K++;
	}
    }
}


$nx=$NZ;
$ny=$K;
$pdl=zeroes($nx,$K);
$k=0;
for ($i=0;$i<$NX;$i++) {
    for ($j=0;$j<$NY;$j++) {
	$val= $pdl_cube->at($i,$j,$nz);
	if ($val>0) {
	    $t=$pdl->slice(",($k)");
	    $t .= $pdl_cube->slice("($i),($j),");
	    $k++
	}
    }
}





$w_min=1e12;
$w_max=-1e12;
$ni=0;
open(DATA,"<$e_file");
while($line=<DATA>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$name[$ni]=$data[1];
	$W[$ni]=$data[0];
	$L1[$ni]=$data[0]-0.5*$width;
	$L2[$ni]=$data[0]+0.5*$width;
	$Lb1[$ni]=$data[0]-2*$width;
	$Lb2[$ni]=$data[0]-1*$width;
	$Lr1[$ni]=$data[0]+1*$width;
	$Lr2[$ni]=$data[0]+2*$width;
	$Lb[$ni]=($Lb1[$ni]+$Lb2[$ni])/2;
	$Lr[$ni]=($Lr1[$ni]+$Lr2[$ni])/2;
	if ($Lb1[$ni]<$w_min) {
	    $w_min=$Lb1[$ni];
	}
	if ($Lr2[$ni]>$w_max) {
	    $w_max=$Lr2[$ni];
	}
	$ni++;
    }
}
close(DATA);
$w_min=$w_min-200;
$w_max=$w_max+200;


$min=1e12;
$max=-1e12;
for ($J=0;$J<$ny;$J++) {
    for ($i=0;$i<$nx;$i++) {
	$w[$i]=$crval+$cdelt*$i;
	$id[$i]=$i;
	$f[$i]=$pdl->at($i,$J);
	if ($min>$f[$i]) {
	    $min=$f[$i];
	}
	if ($max<$f[$i]) {
	    $max=$f[$i];
	}
    }
    $n=$nx;

    pgbegin(0,$dev,1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
#pgenv($w[0],$w[$n-1],$min-5,$max+5,0,0);
    pgenv($w_min,$w_max,-1,1,0,0);
    pgsch(1.6);           # Set character height
    pglabel("Wavelength","Flux","");
    pgsch(2.2);           # Set character height
    pgsci(1);
#pgline($n,\@w,\@f);    
    pgbin($n,\@w,\@f,1);    

    print "$spec_file ";
#    $cdelt=$w[1]-$w[0];
    for ($j=0;$j<$ni;$j++) {
	print "$name[$j] ";
	$Sb=0;
	$Sr=0;
	$nr=0;
	$nb=0;
	$i1_b=($Lb1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	$i2_b=($Lb2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	for ($i=int($i1_b+1);$i<int($i2_b);$i++) {
	    $Sb=$Sb+$f[$i];
#	    print "$i $f[$i] $Sb\n";
	    $nb++;
	}
	$Sb=$Sb/$nb;
	
	$i1_r=($Lr1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	$i2_r=($Lr2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	for ($i=int($i1_r+1);$i<int($i2_r);$i++) {
	    $Sr=$Sr+$f[$i];
#	    print "$i $f[$i] $Sr '$Sb'\n";
	    $nr++;
	}
	$Sr=$Sr/$nr;
	
	$EW=0;
	$k=0;
	$i1=($L1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	$i2=($L2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
	$C=0.5*($Sb+$Sr);
	for ($i=int($i1+1);$i<int($i2);$i++) {
	    $EW=$EW+($f[$i]-$C);
#	    print "$i $f[$i] $Sr $Sb $C $EW\n";
	    $k++;
	}

	print " $EW ";
	pgsls(2);
	pgsci(5);
	pgsfs(2);
	pgpoly(4,[$Lb1[$j],$Lb1[$j],$Lb2[$j],$Lb2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
#    pgline(2,[$Lb2[$j],$Lb2[$j]],[$max*0.5,$max*0.8]);
	pgsci(2);
	pgpoly(4,[$Lr1[$j],$Lr1[$j],$Lr2[$j],$Lr2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
#    pgline(2,[$Lr1[$j],$Lr1[$j]],[$max*0.5,$max*0.8]);
#    pgline(2,[$Lr2[$j],$Lr2[$j]],[$max*0.5,$max*0.8]);
	pgsls(1);
	pgsci(8);
	pgline($k,\@wk,\@CK);    
	pgpoint(1,[$Lb[$j]],[$Sb],2);
	pgpoint(1,[$Lr[$j]],[$Sr],2);
	pgsci(3);
	pgpoly(4,[$L1[$j],$L1[$j],$L2[$j],$L2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
	pgsci(1);
    }
	pgclose;
	pgend;

    print "$J/$ny\n";
    if ($dev eq "/xs") {
	print "Press Enter\n"; <stdin>;
    }
}




exit;
